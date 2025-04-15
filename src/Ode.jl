
module Ode
export daedalus_ode!

using ..Constants
using ..DaedalusStructs

"""
    daedalus_ode!(du, u, p, t)

The ODE system for the DAEDALUS model. This function is intended to be called
    internally from `daedalus`.
"""
function daedalus_ode!(du::Array, u::Array, p::Params, t::Number)
    # du auto-magically takes the dims and type of u (?)
    # each element of the tuple is one of the required params
    # time dependent vaccination
    nu = p.nu * p.switch

    # view the values of each compartment per age group
    # rows represent age groups, epi compartments are columns
    S = @view u[:, iS, :]
    E = @view u[:, iE, :]
    Is = @view u[:, iIs, :]
    Ia = @view u[:, iIa, :]
    H = @view u[:, iH, :]
    R = @view u[:, iR, :]

    # calculate new infections and re-infections
    community_infectious = sum(Is .+ Ia * p.epsilon, dims=2)
    foi = p.beta * p.contacts * community_infectious

    # NOTE: element-wise multiplication
    new_I = S .* foi

    # workplace
    new_I_work = S[i_ECON_GROUPS, :] .* p.cw .*
                 (Is[i_ECON_GROUPS, :] .+ (Ia[i_ECON_GROUPS, :] * p.epsilon))

    new_I[i_ECON_GROUPS, :] += new_I_work

    # views to the change array slice
    dS = @view du[:, 1, :]
    dE = @view du[:, 2, :]
    dIs = @view du[:, 3, :]
    dIa = @view du[:, 4, :]
    dH = @view du[:, 5, :]
    dR = @view du[:, 6, :]
    dD = @view du[:, 7, :]

    # new vaccinations from susceptibles
    new_Svax = @view S[:, 1]
    new_Swane = @view S[:, 2]

    # new vaccinations from recovered
    new_Rvax = @view R[:, 1]
    new_Rwane = @view R[:, 2]

    # calculate change in compartment size and update the change matrix dU
    # note the use of @. for broadcasting, equivalent to .=
    # change in susceptibles
    @. dS = -new_I + (p.rho * R)
    @. dS[:, 1] += (-new_Svax * p.nu + new_Swane * p.psi)
    @. dS[:, 2] += (new_Svax * p.nu - new_Swane * p.psi)

    # change in exposed
    @. dE = new_I - (p.sigma * E)

    # calculate exposed to Is and Ia
    E_sigma = (p.sigma * E) * p.p_sigma
    E_inv_sigma = (p.sigma * E) * (1.0 - p.p_sigma)

    # change in infectious symptomatic
    @. dIs = E_sigma - ((p.gamma_Is .+ p.eta) .* Is)

    # change in infectious asymptomatic
    @. dIa = E_inv_sigma - (p.gamma_Ia * Ia)

    # change in hospitalised
    @. dH = (p.eta .* Is) - ((p.gamma_H + p.omega) .* H)

    # change in recovered
    @. dR = (p.gamma_Ia * Ia) + (p.gamma_Is * Is) +
            (p.gamma_H .* H) - (p.rho * R)

    @. dR[:, 1] += (-new_Rvax * p.nu + new_Rwane * p.psi)
    @. dR[:, 2] += (new_Rvax * p.nu - new_Rwane * p.psi)

    # change in dead
    @. dD = p.omega .* H
end

end


