
export daedalus

using DifferentialEquations
using LinearAlgebra
using StaticArrays

using .Constants
using .Ode
using .Data
using .Events
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
    beta=0.01, # manual beta
    sigma=0.217,
    p_sigma=0.867,
    epsilon=0.58,
    rho=0.003,
    eta::Vector{Float64}=[0.018, 0.082, 0.018, 0.246],
    omega::Vector{Float64}=[0.012, 0.012, 0.012, 0.012],
    gamma_Ia=0.476,
    gamma_Is=0.25,
    gamma_H::Vector{Float64}=[0.034, 0.034, 0.034, 0.034],
    nu=(2.0 / 100.0) / 7.0,
    psi::Float64=1.0 / 270.0,
    npi::Union{Npi,Nothing}=nothing,
    time_end::Float64=100.0,
    increment::Float64=1.0)

    # scale contacts by demography; divide col-wise
    contacts = contacts ./ demography

    eta = [eta; repeat([eta[i_WORKING_AGE]], N_ECON_GROUPS)]
    omega = [omega; repeat([omega[i_WORKING_AGE]], N_ECON_GROUPS)]
    gamma_H = [gamma_H; repeat([gamma_H[i_WORKING_AGE]], N_ECON_GROUPS)]

    # combined parameters into an array; this is not recommended but this cannot be a tuple
    # using a StaticArray for the `contacts` helps cut computation as this is assigned only once(?)
    parameters = Params(contacts, cw, beta, beta, sigma, p_sigma, epsilon,
        rho, eta, omega, omega, gamma_Ia, gamma_Is, gamma_H, nu, psi)

    # prepare the timespan
    timespan = (0.0, time_end)

    # define the ode problem
    ode_problem = ODEProblem(
        daedalus_ode!, initial_state, timespan, parameters
    )

    # check if NPI is passed and define a callback if so
    if (isnothing(npi))
        cb_set = CallbackSet()
        cb_times = []
    else
        fn_effect_on = make_param_changer("beta", *, 0.2)
        fn_effect_off = make_param_reset("beta")

        cb_set = make_time_events(npi, fn_effect_on, fn_effect_off)
        cb_times = get_times(npi)
    end

    # get the solution, ensuring that tstops includes t_vax
    ode_solution = solve(ode_problem, save_everystep=true,
        callback=cb_set,
        tstops=cb_times
    )

    return ode_solution
end
