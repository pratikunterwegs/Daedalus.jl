
module DaedalusStructs

using StaticArrays

export Params, ResponseData, Npi

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
    ResponseData

A struct to hold common response parameters.
"""
struct ResponseData
    time_on::Vector{Float64}
    time_off::Vector{Float64}
end

function ResponseData(; time_on::Vector{Float64}, time_off::Vector{Float64})
    # input checking
    if length(time_on) != length(time_off)
        throw(ArgumentError("time_on and time_off must have the same length."))
    end

    return ResponseData(time_on, time_off)
end

struct Npi
    resparams::ResponseData
    params::NamedTuple
end

function Npi(; params::NamedTuple, time_on::Vector{Float64}, time_off::Vector{Float64})
    resparams = ResponseData(time_on=time_on, time_off=time_off)

    return Npi(resparams, params)
end

end
