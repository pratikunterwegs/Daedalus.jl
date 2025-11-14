
module Ode
export daedalus_ode!

using ..Constants
using ..DaedalusStructs
using ..Helpers

using LinearAlgebra

"""
    daedalus_ode!(du, u, p, t)

The ODE system for the DAEDALUS model. This function is intended to be called
    internally from `daedalus`.
"""
function daedalus_ode!(du::Array, u::Array, p::Params, t::Number)
    # du auto-magically takes the dims and type of u (?)
    # each element of the tuple is one of the required params
    # time dependent vaccination
    nu_eff = p.nu

    # view the values of each compartment per age group
    # rows represent age groups, epi compartments are columns
    U = @view u[1:p.size]
    U = reshape(U, (N_TOTAL_GROUPS, N_COMPARTMENTS, N_VACCINE_STRATA))

    dU = @view du[1:p.size]
    dU = reshape(dU, (N_TOTAL_GROUPS, N_COMPARTMENTS, N_VACCINE_STRATA))

    # using views does not seem to affect performance greatly
    S = @view U[:, iS, :]
    E = @view U[:, iE, :]
    Is = @view U[:, iIs, :]
    Ia = @view U[:, iIa, :]
    H = @view U[:, iH, :]
    R = @view U[:, iR, :]

    # calculate new infections and re-infections
    community_infectious = sum(Is .+ Ia * p.epsilon, dims=2)
    foi = p.beta_now * p.contacts * community_infectious

    # NOTE: element-wise multiplication
    new_I = S .* foi

    # views to the change array slice
    dS = @view dU[:, iS, :]
    dE = @view dU[:, iE, :]
    dIs = @view dU[:, iIs, :]
    dIa = @view dU[:, iIa, :]
    dH = @view dU[:, iH, :]
    dR = @view dU[:, iR, :]
    dD = @view dU[:, iD, :]

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
    @. dS[:, 1] += (-new_Svax * nu_eff + new_Swane * p.psi)
    @. dS[:, 2] += (new_Svax * nu_eff - new_Swane * p.psi)

    # change in exposed
    @. dE = new_I - (p.sigma * E)

    # change in infectious symptomatic
    @. dIs = (p.sigma * E) * p.p_sigma - ((p.gamma_Is .+ p.eta) .* Is)

    # change in infectious asymptomatic
    @. dIa = (p.sigma * E) * (1.0 - p.p_sigma) - (p.gamma_Ia * Ia)

    # change in hospitalised
    @. dH = (p.eta .* Is) - ((p.gamma_H + p.omega_now) .* H)

    # change in recovered
    @. dR = (p.gamma_Ia * Ia) + (p.gamma_Is * Is) +
            (p.gamma_H .* H) - (p.rho * R)

    @. dR[:, 1] += (-new_Rvax * nu_eff + new_Rwane * p.psi)
    @. dR[:, 2] += (new_Rvax * nu_eff - new_Rwane * p.psi)

    # change in dead
    @. dD = p.omega_now .* H
end

end
