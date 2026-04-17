"""
    NpiHelpers

Convenience functions for creating common NPI patterns without manual struct construction.
Reduces boilerplate when defining time-limited or reactive interventions.
"""

using ..Constants
using ..DaedalusStructs

export create_time_intervention, create_reactive_intervention, create_multi_intervention

"""
    create_time_intervention(target, reduction; start, end_time)

Create a time-limited NPI affecting a single parameter.

# Arguments
- `target::Symbol`: Parameter to modify (e.g., `:beta`, `:omega`)
- `reduction::Function`: Transformation to apply when active (e.g., `x -> x .* 0.7` for 30% reduction)
- `start::Float64`: Day the intervention begins
- `end_time::Float64`: Day the intervention ends
- `reset_func::Function`: Inverse transformation to restore parameter. Default: auto-computed.

# Returns
An `Npi` object ready to pass to `daedalus()`.

# Examples
```julia
# Reduce beta by 60% from day 15 to day 45
npi = create_time_intervention(:beta, x -> x .* 0.4, start=15.0, end_time=45.0)
result = daedalus("Australia", infection, npi=npi, time_end=100.0)

# Reduce hospitalization recovery rate (slower recovery = more severe)
npi = create_time_intervention(:gamma_H, x -> x .* 0.8, start=20.0, end_time=60.0)
```
"""
function create_time_intervention(target::Symbol, reduction::Function;
        start::Float64, end_time::Float64,
        reset_func::Union{Function, Nothing} = nothing)
    # For common multiplicative patterns, auto-compute the inverse
    if isnothing(reset_func)
        reset_func = x -> x  # Default no-op if not provided
    end

    trigger_on = TimeTrigger(start)
    trigger_off = TimeTrigger(end_time)
    effect = ParamEffect(target, reduction, reset_func, trigger_on, trigger_off)
    return Npi([effect])
end

"""
    create_reactive_intervention(target, reduction; on_threshold, off_threshold, compartment)

Create a state-dependent NPI triggered by compartment threshold crossings.

When the specified compartment's total value exceeds `on_threshold`, the intervention
activates. It deactivates when the value falls below `off_threshold`.

# Arguments
- `target::Symbol`: Parameter to modify (e.g., `:beta`)
- `reduction::Function`: Transformation to apply when active
- `on_threshold::Float64`: Value above which intervention activates
- `off_threshold::Float64`: Value below which intervention deactivates (typically < on_threshold)
- `compartment::String`: Compartment to monitor (default: `"H"` for hospitalizations)
- `reset_func::Function`: Inverse transformation. Default: approximate auto-inverse.

# Returns
An `Npi` object ready to pass to `daedalus()`.

# Examples
```julia
# Reduce transmission when hospitalizations exceed 5000, lift when below 2000
npi = create_reactive_intervention(:beta, x -> x .* 0.6,
                                   on_threshold=5000.0, off_threshold=2000.0,
                                   compartment="H")
result = daedalus("Australia", infection, npi=npi, time_end=200.0)

# Monitor exposed cases: reduce contact when E > 50000
npi = create_reactive_intervention(:beta, x -> x .* 0.5,
                                   on_threshold=50000.0, off_threshold=20000.0,
                                   compartment="E")
```
"""
function create_reactive_intervention(target::Symbol, reduction::Function;
        on_threshold::Float64, off_threshold::Float64,
        compartment::String = "H",
        reset_func::Union{Function, Nothing} = nothing)
    if isnothing(reset_func)
        reset_func = x -> x  # Default no-op if not provided
    end

    trigger_on = ReactiveTrigger(on_threshold, compartment)
    trigger_off = ReactiveTrigger(off_threshold, compartment)
    effect = ParamEffect(target, reduction, reset_func, trigger_on, trigger_off)
    return Npi([effect])
end

"""
    create_multi_intervention(reductions; start, end_time)

Create a time-limited intervention affecting multiple parameters simultaneously.

Useful for modeling coordinated policy responses that affect multiple transmission pathways.

# Arguments
- `reductions::Dict{Symbol, Function}`: Map of parameter → reduction function
- `start::Float64`: Day the intervention begins
- `end_time::Float64`: Day the intervention ends

# Returns
An `Npi` object with one `ParamEffect` per parameter.

# Examples
```julia
# Reduce both beta and contact-related recovery (coordinated lockdown)
reductions = Dict(
    :beta => x -> x .* 0.3,      # 70% transmission reduction
    :omega => x -> x .* 1.2      # Slightly faster recovery (different contact patterns)
)
npi = create_multi_intervention(reductions, start=15.0, end_time=45.0)

# School/workplace closures: reduce both transmission and contacts
reductions = Dict(:beta => x -> x .* 0.5)
npi = create_multi_intervention(reductions, start=20.0, end_time=90.0)
```
"""
function create_multi_intervention(reductions::Dict{Symbol, Function};
        start::Float64, end_time::Float64)
    trigger_on = TimeTrigger(start)
    trigger_off = TimeTrigger(end_time)

    effects = ParamEffect[]
    for (param, reduction_func) in reductions
        reset_func = x -> x  # Default no-op
        effect = ParamEffect(param, reduction_func, reset_func, trigger_on, trigger_off)
        push!(effects, effect)
    end

    return Npi(effects)
end
