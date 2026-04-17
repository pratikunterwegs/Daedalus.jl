
module DaedalusStructs

using ..Constants

using DiffEqCallbacks

export Params, Npi, Vaccination, Effect, ParamEffect, Trigger, ReactiveTrigger, TimeTrigger

"""
    Params
    
A mutable struct that holds the parameters for the DAEDALUS model.
    It contains the contact matrices, contact weights, and various
    epidemiological parameters.
"""
mutable struct Params
    contacts::Array{Float64, 3}
    cm_temp::Matrix{Float64}
    settings::Int
    ngm::Matrix{Float64}
    demography::Vector{Float64}
    cw::Vector{Float64}
    beta::Float64
    beta_now::Float64  # holds current beta
    sigma::Float64
    p_sigma::Float64
    epsilon::Float64
    rho::Float64
    eta::Vector{Float64}
    omega::Vector{Float64}
    omega_now::Vector{Float64}  # holds current omega
    gamma_Ia::Float64
    gamma_Is::Float64
    gamma_H::Vector{Float64}
    nu::Float64
    psi::Float64
    size_state::Int
end

abstract type Event end

abstract type Trigger end

abstract type Effect end

"""
    ReactiveTrigger

A trigger based on epidemic state (e.g., hospitalization count, Rt).

# Fields
- `name::String`: Compartment name (e.g., "H", "Rt")
- `value::Float64`: Threshold value for triggering

# Constructor
```julia
ReactiveTrigger(value::Float64, name::String)
```
"""
struct ReactiveTrigger <: Trigger
    name::String
    value::Float64

    function ReactiveTrigger(value, name)
        return new(name, value)
    end
end

"""
    TimeTrigger

A trigger based on simulation time.

# Fields
- `name::String`: Always "time"
- `value::Float64`: Time point (days) for triggering

# Constructor
```julia
TimeTrigger(value::Float64, name::String = "time")
```
"""
struct TimeTrigger <: Trigger
    name::String
    value::Float64

    function TimeTrigger(value, name = "time")
        return new(name, value)
    end
end

"""
    ParamEffect

A parameter modification triggered by state or time conditions.

# Fields
- `target::Symbol`: Parameter to modify (e.g., `:beta`, `:omega`)
- `func::Function`: Applied when activated
- `reset_func::Function`: Applied when deactivated
- `trigger_on::Trigger`: Activation condition (`ReactiveTrigger` or `TimeTrigger`)
- `trigger_off::Trigger`: Deactivation condition
- `saved_values::SavedValues`: Internal state tracking
- `ison::Bool`: Current activation state

# Constructor
```julia
ParamEffect(target::Symbol, func::Function, reset_func::Function,
            trigger_on::Trigger, trigger_off::Trigger)
```

# Example
```julia
trigger_on = ReactiveTrigger(5000.0, "H")
trigger_off = ReactiveTrigger(1.0, "Rt")
effect = ParamEffect(:beta, x -> x .* 0.4, x -> x ./ 0.4, trigger_on, trigger_off)
```
"""
mutable struct ParamEffect <: Effect
    target::Symbol
    func::Function
    reset_func::Function
    trigger_on::Trigger
    trigger_off::Trigger
    saved_values::SavedValues
    ison::Bool

    function ParamEffect(target::Symbol, func::Function,
            reset_func::Function,
            trigger_on::Trigger, trigger_off::Trigger)
        sv = SavedValues(Float64, Tuple{Float64, Float64})
        return new(target, func, reset_func, trigger_on, trigger_off, sv, false)
    end
end

"""
    Npi

A container for parameter modifications (effects), each with trigger conditions.

# Fields
- `effects::Vector{Effect}`: Vector of parameter modifications to apply

# Constructor
```julia
Npi(effects::AbstractVector{<:Effect})
```

# Example
```julia
trigger_on = ReactiveTrigger(5000.0, "H")
trigger_off = ReactiveTrigger(1.0, "Rt")
effect = ParamEffect(
    :beta, x -> x .* 0.4, x -> x ./ 0.4, trigger_on, trigger_off)
npi = Npi([effect])
```
"""
mutable struct Npi <: Event
    effects::Vector{Effect}

    function Npi(effects::Vector{<:Effect})
        return new(Vector{Effect}(effects))
    end
end

"""
    Vaccination

An event type that generates vaccination callbacks for the ODE solver.
Modifies `nu` (vaccination rate) and `psi` (waning rate) parameters dynamically.

# Fields
- `start_time::Float64`: Day vaccination campaign begins
- `rate::Float64`: Fraction of population vaccinated per day
- `uptake_limit::Float64`: Maximum coverage fraction (0-1) before campaign stops
- `efficacy::Float64`: Vaccine efficacy (stored for analysis; affects psi in ODE)
- `waning_period::Float64`: Days until vaccine wanes (sets `psi = 1/waning_period`)
- `ison::Bool`: Whether campaign is currently active (internal state)
"""
mutable struct Vaccination <: Event
    start_time::Float64
    rate::Float64
    uptake_limit::Float64
    efficacy::Float64
    waning_period::Float64
    ison::Bool

    function Vaccination(start_time, rate, uptake_limit, efficacy, waning_period)
        return new(start_time, rate, uptake_limit, efficacy, waning_period, false)
    end
end

end
