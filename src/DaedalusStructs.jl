
module DaedalusStructs

using StaticArrays

export Params

mutable struct Params
    contacts::StaticArray
    cw::Vector{Float64}
    beta::Float64
    sigma::Float64
    p_sigma::Float64
    epsilon::Float64
    rho::Float64
    eta::Vector{Float64}
    omega::Vector{Float64}
    gamma_Ia::Float64
    gamma_Is::Float64
    gamma_H::Vector{Float64}
    nu::Float64
    psi::Float64
    t_vax::Float64
    switch::Bool
end
end

