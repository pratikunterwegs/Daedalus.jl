
export daedalus, prepare_shared_data, daedalus_internal, extract_infection_params

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
    daedalus(country, infection; npi=nothing, log_rt=true, time_end=100.0, increment=1.0, n_threads=1)::Union{NamedTuple, Vector{NamedTuple}}

Model the progression of an epidemic with infection parameters and optional interventions.

Solves the DAEDALUS compartmental epidemiological model for a specified country
and infection parameters, with optional non-pharmaceutical interventions.

# Arguments
- `country`: Country to model. Can be:
  - A `String` country name (e.g. `"Australia"`), looked up via `DataLoader.get_country`
  - A [`DataLoader.CountryData`](@ref) struct for custom or pre-fetched data
- `infection`: Infection parameters. Can be:
  - A `String` pathogen name (e.g. `"COVID-19"`), looked up via `DataLoader.get_pathogen`
  - A [`DataLoader.InfectionData`](@ref) struct with all epidemiological parameters
  - A `Vector{DataLoader.InfectionData}` for ensemble runs with multiple pathogens (R0s extracted)
- `npi`: Non-pharmaceutical intervention. Can be:
  - `Npi`: Reactive NPI that responds to epidemic state (hospitalizations, Rt)
  - `TimedNpi`: Time-limited NPI with predefined start/end times
  - `nothing`: No intervention (default)
- `log_rt`: Whether to compute and log effective reproduction number (default: true)
- `time_end`: Simulation end time in days (default: 100.0)
- `increment`: Savepoint interval in days (default: 1.0)
- `n_threads`: Number of threads for multi-run execution (default: 1, ignored for scalar runs)

# Examples

## String pathogen name
```julia
result = daedalus("Australia", "sars-cov-2 delta", time_end=200.0)
```

## Using a pre-fetched CountryData struct
```julia
cd = Daedalus.DataLoader.get_country("Australia")
result = daedalus(cd, "sars-cov-2 delta", time_end=200.0)
```

## InfectionData object with customization
```julia
infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
infection.r0 = 2.5  # Customize R0 if needed
result = daedalus("Australia", infection, time_end=200.0)
```

## Time-limited intervention
```julia
npi = TimedNpi(15.0, 45.0, 0.7, "moderate_lockdown")
result = daedalus("Australia", "sars-cov-2 delta", time_end=200.0, npi=npi)
```

## Multiple pathogens with multi-threading
```julia
infections = [
    Daedalus.DataLoader.get_pathogen("sars-cov-2 delta"),
    Daedalus.DataLoader.get_pathogen("influenza 2009")
]
results = daedalus("Australia", infections, time_end=200.0, n_threads=4)
```

# Returns
- Scalar infection (String or InfectionData): A `NamedTuple` with fields:
  - `sol`: ODE solution object
  - `saves`: Saved values from reactive NPI (or nothing)
  - `npi`: The NPI used (or nothing)
- Vector infection: A `Vector{NamedTuple}`, one per infection, with additional field:
  - `r0`: The r0 value for that run
"""

"""
    extract_infection_params(infection::DataLoader.InfectionData)::NamedTuple

Extract and process epidemiological parameters from an InfectionData object.

Converts base parameters to model-compatible form:
- Expands age-varying parameters (eta, omega, gamma_H) to work with both
  age groups (4) and economic groups (45), matching the 49-component state vector
- Computes working-age indices for proper parameter replication

# Arguments
- `infection::DataLoader.InfectionData`: Infection parameters to extract

# Returns
A `NamedTuple` with fields:
- `r0::Float64`: Basic reproduction number
- `sigma::Float64`: Incubation rate (1/latent period)
- `p_sigma::Float64`: Proportion symptomatic (presymptomatic transmission scaling)
- `epsilon::Float64`: Relative transmissibility of asymptomatic cases
- `rho::Float64`: Rate of infectiousness decay in hospitalized
- `eta::Vector{Float64}`: Hospitalization rate (expanded to 49 groups)
- `omega::Vector{Float64}`: Recovery rate for asymptomatic (expanded to 49 groups)
- `gamma_Ia::Float64`: Recovery rate from asymptomatic infection
- `gamma_Is::Float64`: Recovery rate from symptomatic infection
- `gamma_H::Vector{Float64}`: Recovery rate for hospitalized (expanded to 49 groups)
- `nu::Float64`: Fixed at 0.0 (not user-configurable)
"""
function extract_infection_params(infection::DataLoader.InfectionData)::NamedTuple
    # Extract base parameters
    r0 = infection.r0
    sigma = infection.sigma
    p_sigma = infection.p_sigma
    epsilon = infection.epsilon
    rho = infection.rho
    gamma_Ia = infection.gamma_Ia
    gamma_Is = infection.gamma_Is

    # Expand age-varying parameters to 49 groups (4 age + 45 economic)
    # All economic groups use working-age parameters
    i_WORKING_AGE = 3  # 20-64 age group (third group, 1-indexed)

    eta_expanded = [infection.eta;
                    repeat([infection.eta[i_WORKING_AGE]], Constants.N_ECON_GROUPS)]
    omega_expanded = [fill(0.012, Constants.N_AGE_GROUPS);
                      repeat([0.012], Constants.N_ECON_GROUPS)]
    gamma_H_val = infection.gamma_H_recovery  # Scalar value
    gamma_H_expanded = [fill(gamma_H_val, Constants.N_AGE_GROUPS);
                        repeat([gamma_H_val], Constants.N_ECON_GROUPS)]

    nu = 0.0  # Not user-configurable

    return (
        r0 = r0,
        sigma = sigma,
        p_sigma = p_sigma,
        epsilon = epsilon,
        rho = rho,
        eta = eta_expanded,
        omega = omega_expanded,
        gamma_Ia = gamma_Ia,
        gamma_Is = gamma_Is,
        gamma_H = gamma_H_expanded,
        nu = nu
    )
end

"""
    prepare_shared_data(cd::DataLoader.CountryData, time_end::Float64,
        increment::Float64)::NamedTuple

Prepare country/structure data that is computed once and reused across all parameter sets.

Expensive operations (contact matrix processing, demographics aggregation) are
performed once rather than repeated for each r0 value in ensemble runs.

# Arguments
- `cd::DataLoader.CountryData`: Country demographic and economic data
- `time_end::Float64`: Final simulation time
- `increment::Float64`: Savepoint interval (typically 1.0 for daily output)

# Returns
A `NamedTuple` with fields:
- `init_state`: Initial epidemic state (with appended Rt slot)
- `contacts_unscaled`: Raw contact matrix
- `cw`: Worker contact rates
- `contacts_array`: 3D contact array for ODE
- `cm_temp`: Temporary contact matrix workspace
- `settings`: Number of contact settings
- `demog`: Population vector
- `timespan`: Tuple (0.0, time_end)
- `savepoints`: Range of savepoint times
"""
function prepare_shared_data(
        cd::DataLoader.CountryData,
        time_end::Float64,
        increment::Float64
)::NamedTuple
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
    daedalus_internal(n_runs, shared_data, param_sets, cb_set)::EnsembleSolution

Internal function that orchestrates multi-run ODE solving using EnsembleProblem.

Uses SciMLBase.EnsembleProblem with EnsembleThreads() solver for efficient
parallel execution across multiple parameter sets. The base ODE problem is
created once, then remade for each trajectory with its corresponding parameters.

# Arguments
- `n_runs::Int`: Number of parameter sets / trajectories
- `shared_data::NamedTuple`: Pre-computed country data from `prepare_shared_data()`
- `param_sets`: Array of Params structs, one per trajectory
- `cb_set::CallbackSet`: Callback set (pre-constructed by caller)

# Returns
An `EnsembleSolution` containing solutions for all trajectories at specified savepoints

# Details
The function constructs an `EnsembleProblem` with:
- A single base ODE problem created from the first parameter set
- A `prob_func` that remakes the problem with the i-th parameter set
- `EnsembleThreads()` solver for multi-threaded execution (respects Julia's thread count)
- Savepoints from `shared_data.savepoints`
"""
function daedalus_internal(
        n_runs::Int,
        shared_data::NamedTuple,
        param_sets,
        cb_set::CallbackSet
)::EnsembleSolution
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

# Method 1: String pathogen name - fetch and delegate to InfectionData method
function daedalus(country::Union{String, DataLoader.CountryData}, infection::String;
        npi::Union{Npi, TimedNpi, Nothing} = nothing,
        log_rt::Bool = true,
        time_end::Float64 = 100.0,
        increment::Float64 = 1.0,
        n_threads::Int = 1)
    infection_data = DataLoader.get_pathogen(lowercase(infection))
    return daedalus(country, infection_data; npi = npi, log_rt = log_rt,
        time_end = time_end, increment = increment, n_threads = n_threads)
end

# Method 2: Single InfectionData object
function daedalus(
        country::Union{String, DataLoader.CountryData}, infection::DataLoader.InfectionData;
        npi::Union{Npi, TimedNpi, Nothing} = nothing,
        log_rt::Bool = true,
        time_end::Float64 = 100.0,
        increment::Float64 = 1.0,
        n_threads::Int = 1)

    # Extract infection parameters
    inf_params = extract_infection_params(infection)

    # Resolve country string to CountryData if needed
    cd = isa(country, String) ? DataLoader.get_country(country) : country

    # Prepare shared data once
    shared_data = prepare_shared_data(cd, time_end, increment)
    shared_data.init_state[end] = inf_params.r0 # add R0 manually

    # Get country data for beta calculation
    contacts_unscaled = total_contacts(prepare_contacts(cd; scaled = false))
    cw = worker_contacts(cd)

    # Calculate beta with single r0
    beta = get_beta(
        contacts_unscaled,
        inf_params.r0, inf_params.sigma, inf_params.p_sigma,
        inf_params.epsilon, inf_params.gamma_Ia, inf_params.gamma_Is
    )

    # NGM
    ngm = get_ngm(
        contacts_unscaled,
        beta, inf_params.sigma, inf_params.p_sigma,
        inf_params.epsilon, inf_params.gamma_Ia, inf_params.gamma_Is
    )

    demog = prepare_demog(cd)
    size::Int = N_TOTAL_GROUPS * N_COMPARTMENTS * N_VACCINE_STRATA
    psi::Float64 = 1.0 / 270.0

    # Create parameters struct
    contacts_array = contacts3d(cd)
    settings = get_settings(cd)
    cm_temp::Matrix{Float64} = ones(N_TOTAL_GROUPS, N_TOTAL_GROUPS)
    parameters::Params = Params(
        contacts_array, cm_temp, settings, ngm, demog,
        cw, beta, beta, inf_params.sigma, inf_params.p_sigma,
        inf_params.epsilon, inf_params.rho, inf_params.eta, inf_params.omega,
        inf_params.omega, inf_params.gamma_Ia, inf_params.gamma_Is,
        inf_params.gamma_H, inf_params.nu, psi, size
    )

    # Build callback set
    savepoints = shared_data.savepoints
    if isnothing(npi)
        cb_set = CallbackSet()
        if log_rt
            rt_logger = make_rt_logger(savepoints)
            cb_set = CallbackSet(rt_logger)
        end
    elseif isa(npi, TimedNpi)
        timed_callbacks = make_timed_npi_callbacks(npi)
        if log_rt
            rt_logger = make_rt_logger(savepoints)
            cb_set = CallbackSet(timed_callbacks, rt_logger)
        else
            cb_set = timed_callbacks
        end
    elseif isa(npi, Npi)
        save_events = make_save_events(npi, savepoints)
        events = make_events(npi, savepoints)
        if log_rt
            rt_logger = make_rt_logger(savepoints)
            cb_set = CallbackSet(events, save_events..., rt_logger)
        else
            cb_set = CallbackSet(events, save_events...)
        end
    end

    # Solve for scalar case
    ode_problem = ODEProblem(
        daedalus_ode!, shared_data.init_state, shared_data.timespan, parameters
    )
    ode_solution = solve(ode_problem, callback = cb_set, saveat = savepoints)

    saved_vals = if isnothing(npi)
        nothing
    elseif isa(npi, Npi)
        [eff.saved_values for eff in npi.effects]
    else
        nothing
    end

    return (sol = ode_solution, saves = saved_vals, npi = npi)
end
