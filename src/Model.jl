
export daedalus

using OrdinaryDiffEq  # consider even lighter deps
using DiffEqCallbacks: SavingCallback, SavedValues, PresetTimeCallback
using LinearAlgebra: eigen
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

# Arguments
- `npi`: Non-pharmaceutical intervention. Can be:
  - `Npi`: Reactive NPI that responds to epidemic state (hospitalizations, Rt)
  - `TimedNpi`: Time-limited NPI with predefined start/end times
  - `nothing`: No intervention (default)

# Examples

## No intervention
```julia
result = daedalus(r0=2.5, time_end=200.0)
```

## Single-phase time-limited intervention
```julia
# 30% transmission reduction from day 15 to day 45
npi = TimedNpi(15.0, 45.0, 0.7, "moderate_lockdown")
result = daedalus(r0=2.5, time_end=200.0, npi=npi)
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
result = daedalus(r0=2.5, time_end=200.0, npi=npi)
```

## Reactive (state-dependent) intervention
```julia
# Triggers when hospitalizations reach threshold, deactivates when Rt < 1
npi = Npi(5000.0, (coef=0.4,))
result = daedalus(r0=2.5, time_end=200.0, npi=npi)
```
"""
function daedalus(;
    initial_state=australia_initial_state(australia_demography()),
    contacts=prepare_contacts(),
    cw=worker_contacts(),
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
    npi::Union{Npi,TimedNpi,Nothing}=nothing,
    log_rt=true,
    time_end::Float64=100.0,
    increment::Float64=1.0)

    # calculate beta
    beta = get_beta(
        prepare_contacts(scaled=false),
        r0, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is
    )

    # age varying parameters
    eta = [eta; repeat([eta[i_WORKING_AGE]], N_ECON_GROUPS)]
    omega = [omega; repeat([omega[i_WORKING_AGE]], N_ECON_GROUPS)]
    gamma_H = [gamma_H; repeat([gamma_H[i_WORKING_AGE]], N_ECON_GROUPS)]

    # NGM
    ngm = get_ngm(
        prepare_contacts(scaled=false),
        r0, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is
    )
    demog = SVector{N_TOTAL_GROUPS}(prepare_demog())

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
    initial_state = [initial_state; r0]

    # define the ode problem
    ode_problem = ODEProblem(
        daedalus_ode!, initial_state, timespan, parameters
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

    # get the solution, ensuring that tstops includes t_vax
    ode_solution = solve(ode_problem, callback=cb_set, saveat=savepoints)

    # Handle saved values - only reactive NPIs have saved_values
    saved_vals = if isnothing(npi)
        nothing
    elseif isa(npi, Npi)
        npi.saved_values
    else  # TimedNpi
        nothing
    end

    return (sol=ode_solution, saves=saved_vals, npi=npi)
end
