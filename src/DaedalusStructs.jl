
module DaedalusStructs

using StaticArrays

export Params, TimedResponseData, TimedNpi, ReactiveResponseData, ReactiveNpi

"""
    Params
    
A mutable struct that holds the parameters for the DAEDALUS model.
It contains the contact matrices, contact weights, and various epidemiological parameters.
"""
mutable struct Params
    contacts::StaticArray
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
    TimedResponseData

A struct to hold common reactive response parameters.
"""
struct ReactiveResponseData
    state_on::Vector{Float64}
    state_off::Vector{Float64}
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

"""
    ReactiveResponseData()

Construct `ReactiveResponseData` structs.
"""
function ReactiveResponseData(; state_on::Vector{Float64}, state_off::Vector{Float64})
    # input checking
    if length(state_on) != length(state_off)
        throw(ArgumentError("state_on and state_off must have the same length."))
    end

    return ReactiveResponseData(state_on, state_off)
end

"""
    TimedNpi

A struct to specify a time-bound NPI.
"""
struct TimedNpi
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
    ReactiveNpi

A struct to specify a reactive NPI that responds to a state value.
"""
struct ReactiveNpi
    resparams::ReactiveResponseData
    params::NamedTuple
end

"""
    ReactiveNpi

Function to specify a reactive NPI that responds to a state value.
"""
function ReactiveNpi(; params::NamedTuple, state_on::Vector{Float64},
    state_off::Vector{Float64})
    resparams = ReactiveResponseData(state_on=state_on, state_off=state_off)

    return ReactiveNpi(resparams, params)
end

end
