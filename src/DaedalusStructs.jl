
module DaedalusStructs

using ..Constants

using DiffEqCallbacks

export Params, Npi, Effect, ParamEffect, Trigger, ReactiveTrigger, TimeTrigger

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

struct ReactiveTrigger <: Trigger
    name::String
    value::Float64

    function ReactiveTrigger(value, name)
        return new(name, value)
    end
end

struct TimeTrigger <: Trigger
    name::String
    value::Float64

    function TimeTrigger(value, name = "time")
        return new(name, value)
    end
end

"""
    ParamEffect
"""
mutable struct ParamEffect <: Effect
    target::Symbol
    func::Function
    reset_func::Function
    trigger_on::Trigger
    trigger_off::Trigger
    saved_values::SavedValues
    ison::Bool

    # Full constructor with StateData objects
    function ParamEffect(target::Symbol, func::Function,
            reset_func::Function,
            trigger_on::Trigger, trigger_off::Trigger)
        sv = SavedValues(Float64, Tuple{Float64, Float64})
        return new(target, func, reset_func, trigger_on, trigger_off, sv, false)
    end
end

"""
    Npi

A mutable struct to specify non-pharmaceutical interventions (NPIs). It is a unified container
of parameter effects (both reactive and timed), each specifying how to modify an ODE parameter
and when to activate/deactivate.

# Fields
- `effects::Vector{ParamEffect}`: Parameter modifications to apply.

# Constructors

## Direct container constructor
```julia
Npi(effects::AbstractVector{<:ParamEffect})
```

# Examples

"""
mutable struct Npi <: Event
    effects::Vector{Effect}

    function Npi(effects::Vector{<:Effect})
        return new(Vector{Effect}(effects))
    end
end

end
