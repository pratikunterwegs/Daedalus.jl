
export daedalus

using DifferentialEquations
using LinearAlgebra
using StaticArrays

using .Constants
using .Ode
using .Data
using .Events
using .Helpers
using .DaedalusStructs

"""
    daedalus()

Model the progression of a daedalus epidemic with multiple optional vaccination
strata.
"""
function daedalus(;
    initial_state=australia_initial_state(australia_demography()),
    contacts=prepare_contacts(),
    cw=worker_contacts(),
    demography=prepare_demog(),
    hospital_capacity::Float64=1000.0,
    r0=1.3, # manual beta assumes R0 = 1.3, infectious period = 7 days
    sigma=0.217,
    p_sigma=0.867,
    epsilon=0.58,
    rho=0.003,
    eta::Vector{Float64}=[0.018, 0.082, 0.018, 0.246],
    omega::Vector{Float64}=[0.012, 0.012, 0.012, 0.012],
    gamma_Ia=0.476,
    gamma_Is=0.25,
    gamma_H::Vector{Float64}=[0.034, 0.034, 0.034, 0.034],
    nu=0.0,
    psi::Float64=1.0 / 270.0,
    npi::Union{TimedNpi,ReactiveNpi,Nothing}=nothing,
    time_end::Float64=100.0,
    increment::Float64=1.0)

    # calculate beta
    beta = get_beta(australia_contacts(), r0, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is)

    # age varying parameters
    eta = [eta; repeat([eta[i_WORKING_AGE]], N_ECON_GROUPS)]
    omega = [omega; repeat([omega[i_WORKING_AGE]], N_ECON_GROUPS)]
    gamma_H = [gamma_H; repeat([gamma_H[i_WORKING_AGE]], N_ECON_GROUPS)]

    # NGM
    ngm = get_ngm(
        australia_contacts(), r0, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is
    )
    demog = SVector{N_AGE_GROUPS}(australia_demography())

    size::Int = N_TOTAL_GROUPS * N_COMPARTMENTS * N_VACCINE_STRATA

    # combined parameters into an array; this is not recommended but this cannot be a tuple
    # using a StaticArray for the `contacts` helps cut computation as this is assigned only once(?)
    parameters = Params(contacts, ngm, demog, cw, beta, beta, sigma, p_sigma,
        epsilon, rho, eta, omega, omega, gamma_Ia, gamma_Is, gamma_H, nu, psi,
        size)

    # prepare the timespan and savepoints
    timespan = (0.0, time_end)
    savepoints = 0.0:increment:time_end

    # add Rt compartment at end
    initial_state = reshape(initial_state, length(initial_state))
    initial_state = [initial_state; r0]  # assume tinf = 7.0

    # saving callback for Rt
    saved_values = SavedValues(Float64, Float64)

    idx_susceptible = Constants.get_indices("S")

    # See implementation of saving derivatives in https://discourse.julialang.org/t/wrong-derivatives-saved-with-savingcallback/92471
    savingcb = SavingCallback(
        (u, t, integrator) -> (_u = similar(u); get_du!(_u, integrator); return last(_u)),
        saved_values, saveat=savepoints)

    function affect!(integrator)
        # assumes saved_values in scope --- to be reconsidered
        if (length(saved_values.saveval) > 0)
            Rt = last(saved_values.saveval)
            if (Rt < 1.0)
                # println("time = $(integrator.t); Rt = $Rt")
            end
        end
    end

    pstcb = PresetTimeCallback(
        savepoints, affect!
    )

    # define the ode problem
    ode_problem = ODEProblem(
        daedalus_ode!, initial_state, timespan, parameters
    )

    # check if NPI is passed and define a callback if so
    if (isnothing(npi))
        cb_set = CallbackSet(savingcb, pstcb)
        cb_times = []
    elseif (typeof(npi) == Daedalus.DaedalusStructs.TimedNpi)
        coef = get_coef(npi)

        fn_effect_on = make_param_changer("beta", *, coef)
        fn_effect_off = make_param_reset("beta")

        cb_set = make_events(npi, fn_effect_on, fn_effect_off)
        cb_set = CallbackSet(cb_set, savingcb)
        cb_times = get_times(npi)
    else
        coef = get_coef(npi)
        fn_effect_on = make_param_changer("beta", .*, coef)
        fn_effect_off = make_param_reset("beta")

        cb_set = make_events(npi, fn_effect_on, fn_effect_off)
        cb_set = CallbackSet(cb_set, savingcb)
        cb_times = []
    end

    # get the solution, ensuring that tstops includes t_vax
    ode_solution = solve(ode_problem, saveat=savepoints,
        callback=cb_set,
        tstops=cb_times
    )

    return (sol=ode_solution, saves=saved_values)
end
