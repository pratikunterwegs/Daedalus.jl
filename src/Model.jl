
export daedalus

using DifferentialEquations
using LinearAlgebra
using StaticArrays

using .Constants
using .Ode
using .Data
using .Events

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
    beta=0.01, # manual beta
    sigma=0.217,
    p_sigma=0.867,
    epsilon=0.58,
    rho=0.003,
    eta=[0.018, 0.082, 0.018, 0.246],
    omega=[0.012, 0.012, 0.012, 0.012],
    gamma_Ia=0.476,
    gamma_Is=0.25,
    gamma_H=[0.034, 0.034, 0.034, 0.034],
    nu=(2 / 100) / 7,
    psi::Number=1 / 270,
    t_vax::Number=200.0,
    time_end::Number=300.0,
    threshold::Number=1000.0,
    increment::Number=1.0)

    # scale contacts by demography; divide col-wise
    contacts = contacts ./ demography

    eta = [eta; repeat([eta[i_WORKING_AGE]], N_ECON_GROUPS)]
    omega = [omega; repeat([omega[i_WORKING_AGE]], N_ECON_GROUPS)]
    gamma_H = [gamma_H; repeat([gamma_H[i_WORKING_AGE]], N_ECON_GROUPS)]

    # combined parameters into an array; this is not recommended but this cannot be a tuple
    # using a StaticArray for the `contacts` helps cut computation as this is assigned only once(?)
    switch = false
    parameters = [contacts, cw, beta, sigma, p_sigma, epsilon,
        rho, eta, omega, gamma_Ia, gamma_Is, gamma_H, nu, psi, t_vax,
        switch]

    # prepare the timespan
    timespan = (0.0, time_end)

    # define the ode problem
    ode_problem = ODEProblem(
        daedalus_ode!, initial_state, timespan, parameters
    )

    cond_beta = make_cond(threshold)
    cb_vax = DiscreteCallback(cond_vax, start_vax!)
    cb_npi = ContinuousCallback(cond_beta, reduce_beta!)

    cb_set = CallbackSet(cb_vax, cb_npi)

    # get the solution
    ode_solution = solve(ode_problem, callback=cb_set, tstops=[t_vax])

    return ode_solution
end
