# Vector InfectionData dispatch for daedalus function
# This file extends the daedalus function from Model.jl to handle Vector{InfectionData} input

"""
    daedalus(country, infections::Vector{InfectionData}; npi, log_rt, time_end, increment, n_threads)

Run the daedalus epidemic model with multiple infection parameters in parallel.

This is the vector-InfectionData dispatch of the `daedalus` function. It orchestrates
multi-run ODE solving across multiple InfectionData objects (extracting R0 values)
using shared computation and optional multi-threading.

# Arguments
- `country`: Country to model (String or CountryData struct)
- `infections::Vector{InfectionData}`: Vector of infection parameter objects
- `npi`: Non-pharmaceutical intervention (optional)
- `log_rt`: Whether to log effective reproduction number (default: true)
- `time_end`: Final simulation time in days (default: 100.0)
- `increment`: Savepoint interval in days (default: 1.0)
- `n_threads`: Number of threads for parallel execution (default: 1)

# Returns
A `Vector{NamedTuple}` where each element contains:
- `sol`: The ODE solution for that infection
- `saves`: Saved event values (if applicable)
- `npi`: The NPI specification used
- `r0`: The R0 value for this infection
"""
function daedalus(country::Union{String, DataLoader.CountryData},
        infections::Vector{DataLoader.InfectionData};
        npi::Union{Npi, Nothing} = nothing,
        log_rt::Bool = true,
        time_end::Float64 = 100.0,
        increment::Float64 = 1.0,
        n_threads::Int = 1)

    # Extract R0 values from infection list
    r0_values = [inf.r0 for inf in infections]

    # Use first infection's other parameters
    inf_params = extract_infection_params(infections[1])

    # Resolve country string to CountryData if needed
    cd = isa(country, String) ? DataLoader.get_country(country) : country

    # Prepare shared data once
    shared_data = prepare_shared_data(cd, time_end, increment)

    # Get country data for beta calculation
    contacts_unscaled = total_contacts(prepare_contacts(cd; scaled = false))
    cw = worker_contacts(cd)
    demog = prepare_demog(cd)
    size::Int = N_TOTAL_GROUPS * N_COMPARTMENTS * N_VACCINE_STRATA
    psi::Float64 = 1.0 / 270.0

    # Create parameter sets - one per r0 value
    contacts_array = contacts3d(cd)
    settings = get_settings(cd)
    param_sets = []

    for r0 in r0_values
        beta = get_beta(
            contacts_unscaled,
            r0, inf_params.sigma, inf_params.p_sigma,
            inf_params.epsilon, inf_params.gamma_Ia, inf_params.gamma_Is
        )

        ngm = get_ngm(
            contacts_unscaled,
            beta, inf_params.sigma, inf_params.p_sigma,
            inf_params.epsilon, inf_params.gamma_Ia, inf_params.gamma_Is
        )

        cm_temp::Matrix{Float64} = ones(N_TOTAL_GROUPS, N_TOTAL_GROUPS)
        params = Params(
            contacts_array, cm_temp, settings, ngm, demog,
            cw, beta, beta, inf_params.sigma, inf_params.p_sigma,
            inf_params.epsilon, inf_params.rho, inf_params.eta, inf_params.omega,
            inf_params.omega, inf_params.gamma_Ia, inf_params.gamma_Is,
            inf_params.gamma_H, inf_params.nu, psi, size
        )
        push!(param_sets, params)
    end

    # Build callback set
    savepoints = shared_data.savepoints
    if isnothing(npi)
        cb_set = CallbackSet()
        if log_rt
            rt_logger = make_rt_logger(savepoints)
            cb_set = CallbackSet(rt_logger)
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

    # Solve ensemble
    ensemble_solution = daedalus_internal(
        length(r0_values), shared_data, param_sets, cb_set)

    # Format results: Vector of NamedTuples, one per r0
    results = []
    for i in eachindex(r0_values)
        saved_vals = if isnothing(npi)
            nothing
        elseif isa(npi, Npi)
            [eff.saved_values for eff in npi.effects if isa(eff, ReactiveEffect)]
        else
            nothing
        end

        push!(results, (
            sol = ensemble_solution[i], saves = saved_vals, npi = npi, r0 = r0_values[i]))
    end

    return results
end
