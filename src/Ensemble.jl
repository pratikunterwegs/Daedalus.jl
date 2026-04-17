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
        vaccination::Union{Vaccination, Nothing} = nothing,
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

    # Build unified event vector
    all_events::Vector{DaedalusStructs.Event} = DaedalusStructs.Event[]
    !isnothing(npi) && push!(all_events, npi)
    !isnothing(vaccination) && push!(all_events, vaccination)

    # Build callback set
    savepoints = shared_data.savepoints
    if isempty(all_events)
        cb_set = CallbackSet()
        if log_rt
            rt_logger = make_rt_logger(savepoints)
            cb_set = CallbackSet(rt_logger)
        end
    else
        event_cbs = make_events(all_events, savepoints)
        save_cbs = make_save_events(all_events, savepoints)
        if log_rt
            rt_logger = make_rt_logger(savepoints)
            cb_set = CallbackSet(event_cbs, save_cbs..., rt_logger)
        else
            cb_set = CallbackSet(event_cbs, save_cbs...)
        end
    end

    # Solve ensemble
    ensemble_solution = daedalus_internal(
        length(r0_values), shared_data, param_sets, cb_set)

    # Format results: Vector of DaedalusOutput, one per infection
    results = []
    for i in eachindex(r0_values)
        saved_vals = if isempty(all_events)
            nothing
        else
            [eff.saved_values for ev in all_events if isa(ev, Npi)
             for eff in ev.effects if isa(eff, ParamEffect)]
        end

        push!(results, DaedalusOutput(
            ensemble_solution[i], saved_vals, all_events, cd, infections[i], time_end))
    end

    return results
end
