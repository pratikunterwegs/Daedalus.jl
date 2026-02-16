
module DaedalusStructs

using ..Constants

using DiffEqCallbacks
using StaticArrays

export Params, NpiData, Npi, TimedNpi, StateData, CtStateData, DsStateData,
    get_indices, get_comp_threshold, n_phases, total_duration

"""
    Params
    
A mutable struct that holds the parameters for the DAEDALUS model.
    It contains the contact matrices, contact weights, and various
    epidemiological parameters.
"""
mutable struct Params
    contacts::StaticArray
    ngm::StaticArray
    demography::StaticVector
    cw::StaticVector
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
    size::Int
end

abstract type Event end

abstract type StateData end

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

function get_indices(x::StateData)
    return Constants.get_indices(x.name)
end

function get_comp_threshold(x::StateData)
    return x.value
end

"""
    NpiData

A struct to hold common reactive response parameters.
"""
struct NpiData
    comp_on::StateData
    comp_off::StateData

    function NpiData(value_comp_on::Float64)
        comp_on = CtStateData("H", value_comp_on)
        comp_off = DsStateData("Rt", 1.0)
        return new(comp_on, comp_off)
    end
end

"""
    Npi

A struct to specify an NPI. The idea here is for the NPI to be reactive to
    compartmental values which are logged at discrete time points.
    The struct is mutable so that `saved_values` can be updated.
"""
mutable struct Npi <: Event
    params::NpiData
    coefs::NamedTuple
    saved_values::SavedValues
    ison::Bool

    function Npi(value_on::Float64, coefs::NamedTuple)
        params = NpiData(value_on)
        sv = SavedValues(Float64, Tuple{Bool,Float64,Float64}) # time and two data
        return new(params, coefs, sv, false)
    end
end

"""
    TimedNpi

A struct to specify time-limited non-pharmaceutical interventions (NPIs).
Unlike `Npi`, these interventions are purely time-based and do not respond
to epidemic state. Supports multiple intervention phases with different
intensity levels.

# Fields
- `start_times::Vector{Float64}`: Start times for each intervention phase (days)
- `end_times::Vector{Float64}`: End times for each intervention phase (days)
- `coefs::Vector{Float64}`: Transmission reduction coefficient for each phase.
  Values should be in [0, 1] where 1 = no intervention, 0 = complete transmission block.
- `identifier::String`: Name or description of the intervention strategy

# Constraints
- All vectors must have the same length
- `end_times[i] >= start_times[i]` for all i
- Time intervals must be non-overlapping
- Times must be sorted in ascending order
- Coefficients must be in [0, 1]

# Example
```julia
# Three-phase intervention: moderate → strict → relaxed
timed_npi = TimedNpi(
    [10.0, 30.0, 60.0],  # start times
    [25.0, 55.0, 90.0],  # end times
    [0.7, 0.3, 0.5],     # coefs (0.7 = 30% reduction, 0.3 = 70% reduction)
    "three_phase_lockdown"
)
```
"""
struct TimedNpi <: Event
    start_times::Vector{Float64}
    end_times::Vector{Float64}
    coefs::Vector{Float64}
    identifier::String

    function TimedNpi(
        start_times::Vector{Float64},
        end_times::Vector{Float64},
        coefs::Vector{Float64},
        identifier::String = "custom_timed"
    )
        # Validate vector lengths match
        n_phases = length(start_times)
        if length(end_times) != n_phases
            throw(ArgumentError(
                "end_times must have same length as start_times " *
                "(got $(length(end_times)), expected $n_phases)"
            ))
        end
        if length(coefs) != n_phases
            throw(ArgumentError(
                "coefs must have same length as start_times " *
                "(got $(length(coefs)), expected $n_phases)"
            ))
        end

        # Validate times are non-negative
        if !all(start_times .>= 0.0)
            throw(ArgumentError("start_times must be non-negative"))
        end
        if !all(end_times .>= 0.0)
            throw(ArgumentError("end_times must be non-negative"))
        end

        # Validate end_times >= start_times
        for i in 1:n_phases
            if end_times[i] < start_times[i]
                throw(ArgumentError(
                    "end_times[$i] ($(end_times[i])) must be >= " *
                    "start_times[$i] ($(start_times[i]))"
                ))
            end
        end

        # Validate times are sorted
        if !issorted(start_times)
            throw(ArgumentError("start_times must be in ascending order"))
        end
        if !issorted(end_times)
            throw(ArgumentError("end_times must be in ascending order"))
        end

        # Validate non-overlapping intervals
        for i in 1:(n_phases - 1)
            if end_times[i] >= start_times[i + 1]
                throw(ArgumentError(
                    "Overlapping intervals detected: phase $i ends at " *
                    "$(end_times[i]) but phase $(i+1) starts at " *
                    "$(start_times[i+1])"
                ))
            end
        end

        # Validate coefficients in [0, 1]
        if !all(0.0 .<= coefs .<= 1.0)
            throw(ArgumentError("coefs must be in range [0.0, 1.0]"))
        end

        new(start_times, end_times, coefs, identifier)
    end

    # Convenience constructor for single-phase intervention
    function TimedNpi(
        start_time::Float64,
        end_time::Float64,
        coef::Float64,
        identifier::String = "single_phase"
    )
        TimedNpi([start_time], [end_time], [coef], identifier)
    end
end

"""
    n_phases(npi::TimedNpi)

Return the number of intervention phases in a TimedNpi.

# Example
```julia
npi = TimedNpi([10.0, 30.0], [20.0, 40.0], [0.7, 0.5])
n_phases(npi)  # Returns: 2
```
"""
n_phases(npi::TimedNpi) = length(npi.start_times)

"""
    total_duration(npi::TimedNpi)

Calculate the total duration of all intervention phases combined (not including
gaps between phases).

# Returns
- Total days of active intervention

# Example
```julia
npi = TimedNpi([10.0, 30.0], [20.0, 40.0], [0.7, 0.5])
# Phase 1: 10 days (10-20), Phase 2: 10 days (30-40)
total_duration(npi)  # Returns: 20.0
```
"""
function total_duration(npi::TimedNpi)
    return sum(npi.end_times .- npi.start_times)
end

# Base.show method for pretty printing
function Base.show(io::IO, npi::TimedNpi)
    n = n_phases(npi)
    total_dur = total_duration(npi)

    println(io, "TimedNpi: $(npi.identifier)")
    println(io, "  Number of phases: $n")
    println(io, "  Total duration: $(round(total_dur, digits=1)) days")

    for i in 1:n
        duration = npi.end_times[i] - npi.start_times[i]
        reduction = round((1.0 - npi.coefs[i]) * 100, digits=1)
        println(
            io,
            "  Phase $i: days $(npi.start_times[i])-$(npi.end_times[i]) " *
                "($(round(duration, digits=1))d), $(reduction)% reduction"
        )
    end
end

end
