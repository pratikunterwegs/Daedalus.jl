
module DaedalusStructs

using ..Constants

using DiffEqCallbacks

export Params, Npi, ReactiveEffect, TimedEffect, StateData, CtStateData, DsStateData,
       get_indices

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

abstract type StateData end

abstract type ParamEffect end

struct CtStateData <: StateData
    name::String
    value::Float64
    time_type::String

    function CtStateData(name, value)
        return new(name, value, "continuous")
    end
end

struct DsStateData <: StateData
    name::String
    value::Float64
    time_type::String

    function DsStateData(name, value)
        return new(name, value, "discrete")
    end
end

"""
    ReactiveEffect

A mutable struct specifying how an NPI modifies a single ODE parameter in response to
epidemic state (reactive/state-dependent). Each effect has independent trigger conditions
for activation and deactivation.

# Fields
- `target::Symbol`: The parameter to modify (e.g., `:beta`, `:omega`)
- `func::Function`: A function mapping the original parameter value to the modified value
- `reset_func::Function`: A function mapping the modified parameter value to the
    original value. This is needed when multiple effects overlap, so as to be
    able to return `target` to its pre-modification value.
- `comp_on::StateData`: Trigger condition to activate this effect
- `comp_off::StateData`: Trigger condition to deactivate this effect
- `saved_values::SavedValues`: Historical state values for tracking this effect's triggers
- `is_on::Bool`: Current activation status for this effect

# Constructors

## Full constructor (with StateData objects)
```julia
ReactiveEffect(target::Symbol, func::Function, comp_on::StateData, comp_off::StateData)
```

## Keyword convenience constructor
```julia
ReactiveEffect(target::Symbol, func::Function, reset_func::Function; on::Tuple{String, Float64}, off::Tuple{String, Float64})
```

Example:
```julia
# Reduce beta by 60% when H > 5000, restore when Rt < 1.0
ReactiveEffect(:beta, x -> x .* 0.4, x -> x ./ 0.4; on=("H", 5000.0), off=("Rt", 1.0))

# Reduce gamma_Ia by 20% when D > 100, restore when I < 8000
ReactiveEffect(:gamma_Ia, x -> x .* 0.8, x -> x ./ 0.8; on=("D", 100.0), off=("I", 8000.0))
```
"""
mutable struct ReactiveEffect <: ParamEffect
    target::Symbol
    func::Function
    reset_func::Function
    comp_on::StateData
    comp_off::StateData
    saved_values::SavedValues
    ison::Bool

    # Full constructor with StateData objects
    function ReactiveEffect(target::Symbol, func::Function,
            reset_func::Function,
            comp_on::StateData, comp_off::StateData)
        sv = SavedValues(Float64, Tuple{Float64, Float64})
        return new(target, func, reset_func, comp_on, comp_off, sv, false)
    end
end

# Keyword convenience constructor for ReactiveEffect
function ReactiveEffect(target::Symbol, func::Function, reset_func::Function;
        on::Tuple{String, Float64},
        off::Tuple{String, Float64})
    ReactiveEffect(target, func, reset_func, CtStateData(on...), DsStateData(off...))
end

"""
    get_indices(x::StateData)

Get state vector indices for the compartment in a StateData struct.
Delegates to `Constants.get_indices()`.
"""
function get_indices(x::StateData)
    return Constants.get_indices(x.name)
end

"""
    TimedEffect

A mutable struct specifying how an NPI modifies a single ODE parameter at
specified time points (time-limited / timed intervention).

# Fields
- `target::Symbol`: The parameter to modify (e.g., `:beta`, `:omega`)
- `func::Function`: A function mapping the original parameter value to the modified value
- `reset_func::Function`: A function mapping the modified parameter value back to the original value
- `start_time::Float64`: Time (days) at which to activate the effect
- `end_time::Float64`: Time (days) at which to deactivate the effect

# Constructor
```julia
TimedEffect(target::Symbol, func::Function, reset_func::Function, start_time::Float64, end_time::Float64)
```

# Example
```julia
# Reduce beta by 30% from day 10 to day 40
TimedEffect(:beta, x -> x .* 0.7, x -> x ./ 0.7, 10.0, 40.0)

# Increase omega by 20% from day 15 to day 50
TimedEffect(:omega, x -> x .* 1.2, x -> x ./ 1.2, 15.0, 50.0)
```
"""
mutable struct TimedEffect <: ParamEffect
    target::Symbol
    func::Function
    reset_func::Function
    start_time::Float64
    end_time::Float64

    function TimedEffect(target::Symbol, func::Function, reset_func::Function,
            start_time::Float64, end_time::Float64)
        return new(target, func, reset_func, start_time, end_time)
    end
end

"""
    Npi

A mutable struct to specify non-pharmaceutical interventions (NPIs). It is a unified container
of parameter effects (both reactive and timed), each specifying how to modify an ODE parameter
and when to activate/deactivate.

# Fields
- `effects::Vector{ParamEffect}`: Parameter modifications to apply. May contain both
  `ReactiveEffect` (state-dependent) and `TimedEffect` (time-limited) objects.

# Constructors

## Direct container constructor
```julia
Npi(effects::AbstractVector{<:ParamEffect})
```

# Examples

```julia
# Reactive effects: respond to epidemic state
npi = Npi([
    ReactiveEffect(:beta, x -> x .* 0.4, x -> x ./ 0.4; on=("H", 5000.0), off=("Rt", 1.0)),
    ReactiveEffect(:gamma_Ia, x -> x .* 0.8, x -> x ./ 0.8; on=("D",  100.0), off=("Rt", 1.0)),
])

# Timed effects: activate at specified times
npi = Npi([
    TimedEffect(:beta, x -> x .* 0.7, x -> x ./ 0.7, 10.0, 40.0),
    TimedEffect(:omega, x -> x .* 1.2, x -> x ./ 1.2, 15.0, 50.0),
])

# Mixed: both reactive and timed effects
npi = Npi([
    ReactiveEffect(:beta, x -> x .* 0.6, x -> x ./ 0.6; on=("H", 10000.0), off=("Rt", 1.0)),
    TimedEffect(:omega, x -> x .* 1.1, x -> x ./ 1.1, 20.0, 60.0),
])
```
"""
mutable struct Npi <: Event
    effects::Vector{ParamEffect}

    # Direct constructor (accepts any vector of ParamEffect subtypes)
    function Npi(effects::AbstractVector{<:ParamEffect})
        return new(Vector{ParamEffect}(effects))
    end
end


end
