
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
    t_vax::Float64=200.0,
    npi_times::Vector{Float64}=[10.0, 70.0, 120.0],
    time_end::Float64=300.0,
    hospital_capacity::Float64=1000.0,
    increment::Float64=1.0)

    # scale contacts by demography; divide col-wise
    contacts = contacts ./ demography

    eta = [eta; repeat([eta[i_WORKING_AGE]], N_ECON_GROUPS)]
    omega = [omega; repeat([omega[i_WORKING_AGE]], N_ECON_GROUPS)]
    gamma_H = [gamma_H; repeat([gamma_H[i_WORKING_AGE]], N_ECON_GROUPS)]

    # combined parameters into an array; this is not recommended but this cannot be a tuple
    # using a StaticArray for the `contacts` helps cut computation as this is assigned only once(?)
    switch::Bool = false
    parameters = Params(contacts, cw, beta, sigma, p_sigma, epsilon,
        rho, eta, omega, gamma_Ia, gamma_Is, gamma_H, nu, psi, t_vax,
        switch)

    # prepare the timespan
    timespan = (0.0, time_end)

    # define the ode problem
    ode_problem = ODEProblem(
        daedalus_ode!, initial_state, timespan, parameters
    )

    # make callbacks and callbackset
    condition_hosp_capacity = make_state_condition(hospital_capacity, iH)
    condition_vax_time = make_time_condition(t_vax)
    cb_vax = DiscreteCallback(condition_vax_time, start_vax!)
    cb_npi = DiscreteCallback(make_time_condition(npi_times), reduce_beta!)

    cb_set = CallbackSet(cb_vax, cb_npi)

    # get the solution, ensuring that tstops includes t_vax
    ode_solution = solve(ode_problem, save_everystep=true,
        callback=cb_set, tstops=[t_vax; npi_times])

    return ode_solution
end
