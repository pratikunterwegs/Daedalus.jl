
export daedalus, prepare_shared_data, daedalus_internal

using OrdinaryDiffEq  # consider even lighter deps
using DiffEqCallbacks: SavingCallback, SavedValues, PresetTimeCallback

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
result = daedalus("Australia", 2.5, time_end=200.0)
```

## Using a pre-fetched CountryData struct
```julia
cd = Daedalus.DataLoader.get_country("Australia")
result = daedalus(cd, 2.5, time_end=200.0)
```

## Single-phase time-limited intervention
```julia
# 30% transmission reduction from day 15 to day 45
npi = TimedNpi(15.0, 45.0, 0.7, "moderate_lockdown")
result = daedalus("Australia", 2.5, time_end=200.0, npi=npi)
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
result = daedalus("Australia", 2.5, time_end=200.0, npi=npi)
```

## Reactive (state-dependent) intervention
```julia
# Triggers when hospitalizations reach threshold, deactivates when Rt < 1
npi = Npi(5000.0, (coef=0.4,))
result = daedalus("Australia", 2.5, time_end=200.0, npi=npi)
```

## Multiple R0 values with multi-threading
```julia
# Run model with multiple R0 values in parallel
results = daedalus("Australia", [1.0, 2.0, 3.0], time_end=200.0)
# Returns: Vector of results, one per R0 value
```
"""

"""
    prepare_shared_data(cd::DataLoader.CountryData, time_end::Float64, 
        increment::Float64)

Prepare shared data that is computed once and reused across all parameter sets.
Returns a NamedTuple containing initial state, contact matrices, demographics,
timespan, and callbacks.
"""
function prepare_shared_data(
        cd::DataLoader.CountryData,
        time_end::Float64,
        increment::Float64
)
    # All of these are expensive and independent of infection parameters
    init_state = initial_state(cd)

    # add Rt compartment at end; this holds zero as r0 may differ across runs
    init_state = reshape(init_state, length(init_state))
    init_state = [init_state; 0.0]

    contacts_unscaled = total_contacts(prepare_contacts(cd; scaled = false))
    cw = worker_contacts(cd)
    contacts_array = contacts3d(cd)
    settings = get_settings(cd)
    demog = prepare_demog(cd)
    cm_temp::Matrix{Float64} = ones(N_TOTAL_GROUPS, N_TOTAL_GROUPS)

    # Timespan setup
    timespan = (0.0, time_end)
    savepoints = 0.0:increment:time_end

    return (
        init_state = init_state,
        contacts_unscaled = contacts_unscaled,
        cw = cw,
        contacts_array = contacts_array,
        cm_temp = cm_temp,
        settings = settings,
        demog = demog,
        timespan = timespan,
        savepoints = savepoints
    )
end

"""
    daedalus_internal(n_runs, shared_data, param_sets, cb_set)

Internal function that orchestrates single or multi-run ODE solving.
Handles both scalar r0 and vector r0 cases.
Uses efficient ODE problem reuse via remake() and optional multi-threading.

# Arguments
- `n_runs`: Number of parameter sets
- `shared_data`: Pre-computed country/structure data from prepare_shared_data()
- `param_sets`: Collection of parameter sets.
- `cb_set`: Callback set (pre-constructed by caller)
"""
function daedalus_internal(
        n_runs::Int,
        shared_data::NamedTuple,
        param_sets,
        cb_set::CallbackSet
)
    iRt = Daedalus.Constants.get_indices("Rt")

    # create base problem with initial state
    # Initial state for first run
    init_state_first = shared_data.init_state
    # init_state_first[iRt] = param_sets[1].r0
    params_first = param_sets[1]

    base_problem = ODEProblem(
        daedalus_ode!, init_state_first, shared_data.timespan, params_first
    )

    prob_func = let p = param_sets[1], param_sets = param_sets
        (prob, i, repeat) -> begin
            SciMLBase.remake(base_problem, p = param_sets[i])
        end
    end

    ensemble_prob = EnsembleProblem(base_problem, prob_func = prob_func)

    # NOTE: this does not save data at specific points
    solution = solve(
        ensemble_prob, EnsembleThreads(),
        callback = cb_set,
        trajectories = n_runs,
        saveat = shared_data.savepoints
    )

    return solution
end

function daedalus(country::Union{String, DataLoader.CountryData}, r0::Float64;
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
        increment::Float64 = 1.0)

    # resolve country string to CountryData if needed
    cd = isa(country, String) ? DataLoader.get_country(country) : country

    # get all values from country data
    init_state = initial_state(cd)
    contacts_unscaled = total_contacts(prepare_contacts(cd; scaled = false))
    cw = worker_contacts(cd)

    # calculate beta
    beta = get_beta(
        contacts_unscaled,
        r0, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is
    )

    # age varying parameters
    eta = [eta; repeat([eta[i_WORKING_AGE]], N_ECON_GROUPS)]
    omega = [omega; repeat([omega[i_WORKING_AGE]], N_ECON_GROUPS)]
    gamma_H = [gamma_H; repeat([gamma_H[i_WORKING_AGE]], N_ECON_GROUPS)]

    # NGM
    ngm = get_ngm(
        contacts_unscaled,
        beta, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is
    )
    demog = prepare_demog(cd)

    size::Int = N_TOTAL_GROUPS * N_COMPARTMENTS * N_VACCINE_STRATA

    # combined parameters into an array
    contacts_array = contacts3d(cd)
    settings = get_settings(cd)
    cm_temp::Matrix{Float64} = ones(N_TOTAL_GROUPS, N_TOTAL_GROUPS)
    parameters::Params = Params(contacts_array, cm_temp, settings, ngm, demog,
        cw, beta, beta, sigma, p_sigma, epsilon, rho, eta, omega, omega,
        gamma_Ia, gamma_Is, gamma_H, nu, psi,
        size)

    # prepare the timespan and savepoints
    timespan = (0.0, time_end)
    savepoints = 0.0:increment:time_end

    # add Rt compartment at end
    init_state = reshape(init_state, length(init_state))
    init_state = [init_state; r0]

    # define the ode problem
    ode_problem = ODEProblem(
        daedalus_ode!, init_state, timespan, parameters
    )

    # check if NPI is passed and define callbacks accordingly
    if isnothing(npi)
        # No intervention
        cb_set = CallbackSet()
        if log_rt
            rt_logger = make_rt_logger(savepoints)
            cb_set = CallbackSet(rt_logger)
        end
    elseif isa(npi, TimedNpi)
        # Time-limited NPI - purely time-based triggers
        timed_callbacks = make_timed_npi_callbacks(npi)

        if log_rt
            rt_logger = make_rt_logger(savepoints)
            cb_set = CallbackSet(timed_callbacks, rt_logger)
        else
            cb_set = timed_callbacks
        end
    elseif isa(npi, Npi)
        # Reactive NPI - state-dependent triggers
        coef = get_coef(npi)
        fn_effect_on = make_param_changer("beta", .*, coef)
        fn_effect_off = make_param_reset("beta")

        save_events = make_save_events(npi, savepoints)
        events = make_events(npi, fn_effect_on, fn_effect_off, savepoints)

        if log_rt
            rt_logger = make_rt_logger(savepoints)
            cb_set = CallbackSet(events, save_events, rt_logger)
        else
            cb_set = CallbackSet(events, save_events)
        end
    end

    # get the solution
    ode_solution = solve(ode_problem, callback = cb_set, saveat = savepoints)

    # Handle saved values - only reactive NPIs have saved_values
    saved_vals = if isnothing(npi)
        nothing
    elseif isa(npi, Npi)
        npi.saved_values
    else  # TimedNpi
        nothing
    end

    return (sol = ode_solution, saves = saved_vals, npi = npi)
end
