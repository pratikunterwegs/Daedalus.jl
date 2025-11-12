
module DaedalusStructs

using StaticArrays

export Params, TimedResponseData, TimedNpi, ReactiveResponseData, ReactiveNpi,
    Npi

"""
    Params
    
A mutable struct that holds the parameters for the DAEDALUS model.
It contains the contact matrices, contact weights, and various epidemiological parameters.
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

"""
    TimedResponseData

A struct to hold common timed response parameters.
"""
struct TimedResponseData
    time_on::Vector{Float64}
    time_off::Vector{Float64}
end

"""
    TimedResponseData()

Construct `TimedResponseData` structs.
"""
function TimedResponseData(; time_on::Vector{Float64}, time_off::Vector{Float64})
    # input checking
    if length(time_on) != length(time_off)
        throw(ArgumentError("time_on and time_off must have the same length."))
    end

    return TimedResponseData(time_on, time_off)
end

abstract type Npi end

"""
    TimedNpi

A struct to specify a time-bound NPI.
"""
struct TimedNpi <: Npi
    resparams::TimedResponseData
    params::NamedTuple
end

"""
    TimedNpi

Function to specify a time-bound NPI.
"""
function TimedNpi(; params::NamedTuple, time_on::Vector{Float64}, time_off::Vector{Float64})
    resparams = TimedResponseData(time_on=time_on, time_off=time_off)

    return TimedNpi(resparams, params)
end

"""
    ReactiveResponseData

A struct to hold common reactive response parameters.
"""
struct ReactiveResponseData
    id_state_on::Union{Vector{Int},UnitRange{Int}}
    id_state_off::Union{Vector{Int},UnitRange{Int}}
    value_state_on::Float64
    value_state_off::Float64
end

"""
    ReactiveResponseData()

Construct `ReactiveResponseData` structs.
"""
function ReactiveResponseData(;
    id_state_on::Union{Vector{Int},UnitRange{Int}},
    id_state_off::Union{Vector{Int},UnitRange{Int}},
    value_state_on::Float64,
    value_state_off::Float64)

    # TODO: input checking

    return ReactiveResponseData(
        id_state_on, id_state_off, value_state_on, value_state_off
    )
end

"""
    ReactiveNpi

A struct to specify a reactive NPI that responds to a state value.
"""
struct ReactiveNpi <: Npi
    resparams::ReactiveResponseData
    params::NamedTuple
end

"""
    ReactiveNpi()

Function to specify a reactive NPI that responds to a state value.
"""
function ReactiveNpi(; params::NamedTuple, id_state_on::Union{Vector{Int},UnitRange{Int}},
    id_state_off::Union{Vector{Int},UnitRange{Int}},
    value_state_on::Float64,
    value_state_off::Float64)

    resparams = ReactiveResponseData(
        id_state_on=id_state_on, id_state_off=id_state_off,
        value_state_on=value_state_on, value_state_off=value_state_off
    )

    return ReactiveNpi(resparams, params)
end

end
