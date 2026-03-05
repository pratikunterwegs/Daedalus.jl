
export daedalus

using OrdinaryDiffEq  # consider even lighter deps
using DiffEqCallbacks: SavingCallback, SavedValues, PresetTimeCallback
using LinearAlgebra: eigen
using StaticArrays
using Base.Threads

using .Constants
using .Ode
using .Data
using .DataLoader
using .Events
using .Helpers
using .DaedalusStructs

"""
    daedalus()

Model the progression of a daedalus epidemic with multiple optional vaccination
strata.

# Arguments
- `country`: Country to model. Can be:
  - A `String` country name (e.g. `"Australia"`), looked up via `DataLoader.get_country`
  - A [`DataLoader.CountryData`](@ref) struct for custom or pre-fetched data
- `npi`: Non-pharmaceutical intervention. Can be:
  - `Npi`: Reactive NPI that responds to epidemic state (hospitalizations, Rt)
  - `TimedNpi`: Time-limited NPI with predefined start/end times
  - `nothing`: No intervention (default)

# Examples

## No intervention
```julia
result = daedalus(country="Australia", r0=2.5, time_end=200.0)
```

## Using a pre-fetched CountryData struct
```julia
cd = Daedalus.DataLoader.get_country("Australia")
result = daedalus(country=cd, r0=2.5, time_end=200.0)
```

## Single-phase time-limited intervention
```julia
# 30% transmission reduction from day 15 to day 45
npi = TimedNpi(15.0, 45.0, 0.7, "moderate_lockdown")
result = daedalus(country="Australia", r0=2.5, time_end=200.0, npi=npi)
```

## Multi-phase time-limited intervention
```julia
# Three phases: moderate → strict → relaxed
npi = TimedNpi(
    [10.0, 30.0, 60.0],  # start times
    [25.0, 55.0, 90.0],  # end times
    [0.7, 0.3, 0.5],     # coefficients
    "three_phase_strategy"
)
result = daedalus(country="Australia", r0=2.5, time_end=200.0, npi=npi)
```

## Reactive (state-dependent) intervention
```julia
# Triggers when hospitalizations reach threshold, deactivates when Rt < 1
npi = Npi(5000.0, (coef=0.4,))
result = daedalus(country="Australia", r0=2.5, time_end=200.0, npi=npi)
```

## Multiple R0 values with multi-threading
```julia
# Run model with multiple R0 values in parallel
results = daedalus(country="Australia", r0=[1.0, 2.0, 3.0], time_end=200.0, n_threads=4)
# Returns: Vector of results, one per R0 value
```
"""

"""
    prepare_shared_data(cd::DataLoader.CountryData, time_end::Float64, increment::Float64, log_rt::Bool)

Prepare shared data that is computed once and reused across all parameter sets.
Returns a NamedTuple containing initial state, contact matrices, demographics, timespan, and callbacks.
"""
function prepare_shared_data(
    cd::DataLoader.CountryData,
    time_end::Float64,
    increment::Float64,
    log_rt::Bool
)
    # All of these are expensive and independent of infection parameters
    init_state = initial_state(cd)
    contacts_unscaled = total_contacts(prepare_contacts(cd; scaled = false))
    cw = worker_contacts(cd)
    contacts_array = contacts3d(cd)
    settings = get_settings(cd)
    demog = SVector{N_TOTAL_GROUPS}(prepare_demog(cd))

    # Timespan setup
    timespan = (0.0, time_end)
    savepoints = 0.0:increment:time_end

    return (
        init_state = init_state,
        contacts_unscaled = contacts_unscaled,
        cw = cw,
        contacts_array = contacts_array,
        settings = settings,
        demog = demog,
        timespan = timespan,
        savepoints = savepoints,
    )
end

"""
    daedalus_internal(cd, r0_input, shared_data, sigma, p_sigma, epsilon, rho, eta, omega, gamma_Ia, gamma_Is, gamma_H, nu, psi, cb_set, n_threads)

Internal function that orchestrates single or multi-run ODE solving.
Handles both scalar r0 and vector r0 cases.
Uses efficient ODE problem reuse via remake() and optional multi-threading.

# Arguments
- `r0_input`: Scalar Float64 or Vector{Float64} of r0 values
- `shared_data`: Pre-computed country/structure data from prepare_shared_data()
- `cb_set`: Callback set (pre-constructed by caller)
- `n_threads`: Number of threads for parallel execution
"""
function daedalus_internal(
    cd::DataLoader.CountryData,
    r0_input::Union{Float64, Vector{Float64}},
    shared_data::NamedTuple,
    sigma, p_sigma, epsilon, rho, eta, omega, gamma_Ia, gamma_Is, gamma_H, nu, psi,
    cb_set::CallbackSet,
    n_threads::Int
)
    # Normalize r0_input to vector for uniform handling
    r0_vec = isa(r0_input, Vector) ? r0_input : [r0_input]
    n_runs = length(r0_vec)
    is_scalar_input = !isa(r0_input, Vector)

    results = Vector{NamedTuple}(undef, n_runs)

    # Create base problem for first r0, then reuse structure via remake()
    r0_first = r0_vec[1]

    # Compute derived parameters for first run
    beta_first = get_beta(
        shared_data.contacts_unscaled,
        r0_first, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is
    )
    ngm_first = get_ngm(
        shared_data.contacts_unscaled,
        r0_first, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is
    )

    # Expand age-varying parameters
    eta_exp = [eta; repeat([eta[i_WORKING_AGE]], N_ECON_GROUPS)]
    omega_exp = [omega; repeat([omega[i_WORKING_AGE]], N_ECON_GROUPS)]
    gamma_H_exp = [gamma_H; repeat([gamma_H[i_WORKING_AGE]], N_ECON_GROUPS)]

    size::Int = N_TOTAL_GROUPS * N_COMPARTMENTS * N_VACCINE_STRATA
    cm_temp::Matrix{Float64} = ones(N_TOTAL_GROUPS, N_TOTAL_GROUPS)

    params_first = Params(
        shared_data.contacts_array, cm_temp, shared_data.settings, ngm_first,
        shared_data.demog, shared_data.cw, beta_first, beta_first, sigma, p_sigma,
        epsilon, rho, eta_exp, omega_exp, omega_exp, gamma_Ia, gamma_Is, gamma_H_exp,
        nu, psi, size
    )

    # Initial state for first run (reshape to 1D for concatenation with scalar r0)
    init_state_vec = reshape(shared_data.init_state, length(shared_data.init_state))
    init_state_aug = [init_state_vec; r0_first]

    # Create base ODE problem
    base_problem = ODEProblem(
        daedalus_ode!, init_state_aug, shared_data.timespan, params_first
    )

    # Function to solve a single run (for reuse in both serial and threaded execution)
    function solve_single_run(i::Int)
        r0_i = r0_vec[i]

        # Compute derived parameters for this run
        beta_i = get_beta(
            shared_data.contacts_unscaled,
            r0_i, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is
        )
        ngm_i = get_ngm(
            shared_data.contacts_unscaled,
            r0_i, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is
        )

        # Create params for this run
        params_i = Params(
            shared_data.contacts_array, cm_temp, shared_data.settings, ngm_i,
            shared_data.demog, shared_data.cw, beta_i, beta_i, sigma, p_sigma,
            epsilon, rho, eta_exp, omega_exp, omega_exp, gamma_Ia, gamma_Is,
            gamma_H_exp, nu, psi, size
        )

        # Augment initial state with this run's r0 (reshape to 1D first)
        init_state_i = [init_state_vec; r0_i]

        # Create ODE problem via remake (efficient reuse of base structure)
        problem_i = SciMLBase.remake(base_problem; u0 = init_state_i, p = params_i)

        # Solve
        ode_solution = solve(problem_i, callback = cb_set, saveat = shared_data.savepoints)

        return (sol = ode_solution, saves = nothing, npi = nothing, r0 = r0_i)
    end

    # Execute solves (serial or multi-threaded)
    if n_threads > 1
        # Multi-threaded execution
        Threads.@threads for i in 1:n_runs
            results[i] = solve_single_run(i)
        end
    else
        # Serial execution
        for i in 1:n_runs
            results[i] = solve_single_run(i)
        end
    end

    # Return single result if scalar input, otherwise return vector
    if is_scalar_input
        return results[1]
    else
        return results
    end
end

function daedalus(;
        country::Union{String, DataLoader.CountryData} = "Canada",
        r0::Union{Float64, Vector{Float64}} = 1.3, # manual beta assumes R0 = 1.3, infectious period = 7 days
        sigma = 0.217,
        p_sigma = 0.867,
        epsilon = 0.58,
        rho = 0.003,
        eta::Vector{Float64} = [0.018, 0.082, 0.018, 0.246],
        omega::Vector{Float64} = [0.012, 0.012, 0.012, 0.012],
        gamma_Ia = 0.476,
        gamma_Is = 0.25,
        gamma_H::Vector{Float64} = [0.034, 0.034, 0.034, 0.034],
        nu = 0.0,
        psi::Float64 = 1.0 / 270.0,
        npi::Union{Npi, TimedNpi, Nothing} = nothing,
        log_rt = true,
        time_end::Float64 = 100.0,
        increment::Float64 = 1.0,
        n_threads::Int = 1)

    # resolve country string to CountryData if needed
    cd = isa(country, String) ? DataLoader.get_country(country) : country

    # Prepare shared data (done once for both single and multi-run)
    shared_data = prepare_shared_data(cd, time_end, increment, log_rt)

    # Build callback set (shared across all runs)
    if isnothing(npi)
        if log_rt
            rt_logger = make_rt_logger(shared_data.savepoints)
            cb_set = CallbackSet(rt_logger)
        else
            cb_set = CallbackSet()
        end
    elseif isa(npi, TimedNpi)
        timed_callbacks = make_timed_npi_callbacks(npi)
        if log_rt
            rt_logger = make_rt_logger(shared_data.savepoints)
            cb_set = CallbackSet(timed_callbacks, rt_logger)
        else
            cb_set = timed_callbacks
        end
    elseif isa(npi, Npi)
        coef = get_coef(npi)
        fn_effect_on = make_param_changer("beta", .*, coef)
        fn_effect_off = make_param_reset("beta")
        save_events = make_save_events(npi, shared_data.savepoints)
        events = make_events(npi, fn_effect_on, fn_effect_off, shared_data.savepoints)
        if log_rt
            rt_logger = make_rt_logger(shared_data.savepoints)
            cb_set = CallbackSet(events, save_events, rt_logger)
        else
            cb_set = CallbackSet(events, save_events)
        end
    end

    # Delegate to internal solver for both single and multi-run cases
    if isa(r0, Vector)
        # Multi-run mode
        results = daedalus_internal(
            cd, r0, shared_data,
            sigma, p_sigma, epsilon, rho, eta, omega, gamma_Ia, gamma_Is, gamma_H, nu, psi,
            cb_set, n_threads
        )
        # For multi-run, post-process results to add npi and saved_values
        for i in eachindex(results)
            saved_vals = if isnothing(npi)
                nothing
            elseif isa(npi, Npi)
                npi.saved_values
            else
                nothing
            end
            results[i] = (sol = results[i].sol, saves = saved_vals, npi = npi, r0 = results[i].r0)
        end
        return results
    else
        # Single-run mode - call internal and post-process result
        result = daedalus_internal(
            cd, r0, shared_data,
            sigma, p_sigma, epsilon, rho, eta, omega, gamma_Ia, gamma_Is, gamma_H, nu, psi,
            cb_set, n_threads
        )
        # For single-run, extract saved_values from NPI if reactive
        saved_vals = if isnothing(npi)
            nothing
        elseif isa(npi, Npi)
            npi.saved_values
        else
            nothing
        end
        return (sol = result.sol, saves = saved_vals, npi = npi)
    end
end

