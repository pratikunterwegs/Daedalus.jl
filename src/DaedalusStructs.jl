
module DaedalusStructs

using ..Constants

using DiffEqCallbacks
using StaticArrays

export Params, NpiData, Npi, StateData, CtStateData, DsStateData,
    get_indices, get_comp_threshold

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

end
