# Vector R0 dispatch for daedalus function
# This file extends the daedalus function from Model.jl to handle Vector{Float64} r0 input

"""
    daedalus(; country, r0::Vector{Float64}, ...)

Run the daedalus epidemic model with multiple R0 values in parallel.

This is the vector-R0 dispatch of the `daedalus` function. It orchestrates
multi-run ODE solving across multiple R0 values using shared computation
and optional multi-threading.

# Arguments
- `country`: Country to model (String or CountryData struct)
- `r0::Vector{Float64}`: Vector of reproduction numbers to simulate
- Other arguments: Same as scalar daedalus() function

# Returns
A `Vector{NamedTuple}` where each element contains:
- `sol`: The ODE solution for that R0 value
- `saves`: Saved event values (if applicable)
- `npi`: The NPI specification used
- `r0`: The R0 value for this run
"""
function daedalus(country::Union{String, DataLoader.CountryData},
        r0::Vector{Float64};
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

    n_runs = length(r0)

    # Prepare shared data (done once for all runs)
    shared_data = prepare_shared_data(cd, time_end, increment)

    # prepare betas and NGMs for each R0 value
    beta = get_beta(
        shared_data.contacts_unscaled, r0, sigma, p_sigma,
        epsilon, gamma_Ia, gamma_Is
    )

    ngm = get_ngm(
        shared_data.contacts_unscaled, beta, sigma, p_sigma,
        epsilon, gamma_Ia, gamma_Is
    )

    # Ensure beta and ngm are vectors
    if isa(beta, Float64)
        beta = [beta]
    end
    if !isa(ngm, Vector)
        ngm = [ngm]
    end

    # Expand age-varying parameters
    eta_exp = [eta; repeat([eta[i_WORKING_AGE]], N_ECON_GROUPS)]
    omega_exp = [omega; repeat([omega[i_WORKING_AGE]], N_ECON_GROUPS)]
    gamma_H_exp = [gamma_H; repeat([gamma_H[i_WORKING_AGE]], N_ECON_GROUPS)]

    size_state::Int = N_TOTAL_GROUPS * N_COMPARTMENTS * N_VACCINE_STRATA

    # Create parameter sets for each R0 value
    param_sets = [Params(
                      shared_data.contacts_array, shared_data.cm_temp,
                      shared_data.settings,
                      ngm[i], shared_data.demog,
                      shared_data.cw, beta[i], beta[i],
                      sigma, p_sigma, epsilon, rho, eta_exp,
                      omega_exp, omega_exp, gamma_Ia, gamma_Is,
                      gamma_H_exp, nu, psi, size_state
                  ) for i in 1:n_runs]

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

    # Solve ensemble problem
    solutions = daedalus_internal(
        n_runs, shared_data, param_sets, cb_set
    )

    # return vector of results compatible with current output handlers
    # NOTE: use EnsembleResult handling from SciMLBase later
    results = [(sol = solutions[i], saves = _get_saved_vals_ensemble(npi),
                   npi = npi, r0 = r0[i])
               for i in 1:n_runs]

    return results
end

# Helper function to extract saved values from reactive NPIs
function _get_saved_vals_ensemble(npi::Union{Npi, TimedNpi, Nothing})
    if isnothing(npi)
        return nothing
    elseif isa(npi, Npi)
        return npi.saved_values
    else  # TimedNpi
        return nothing
    end
end
