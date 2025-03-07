
module Ode
export daedalus_ode!

using ..Constants

"""
    daedalus_ode!(du, u, p, t)

The ODE system for the DAEDALUS model. This function is intended to be called
    internally from `daedalus`.
"""
function daedalus_ode!(du, u, p, t)
    # du auto-magically takes the type of u (?)
    # each element of the tuple is one of the required params
    contacts, cw, beta, sigma, p_sigma, epsilon,
    rho, eta, omega, gamma_Ia, gamma_Is, gamma_H, nu, psi,
    t_vax = p

    # time dependent vaccination
    nu = t > t_vax ? nu : 0.0

    # view the values of each compartment per age group
    # rows represent age groups, epi compartments are columns
    S = @view u[:, iS, :]
    E = @view u[:, iE, :]
    Is = @view u[:, iIs, :]
    Ia = @view u[:, iIa, :]
    H = @view u[:, iH, :]
    R = @view u[:, iR, :]

    # calculate new infections and re-infections
    community_infectious = sum(Is .+ Ia * epsilon, dims=2)
    foi = beta * contacts * community_infectious

    # NOTE: element-wise multiplication
    new_I = S .* foi

    # workplace
    new_I_work = S[i_ECON_GROUPS, :] .* cw .*
                 (Is[i_ECON_GROUPS, :] .+ (Ia[i_ECON_GROUPS, :] * epsilon))

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
    @. dS = -new_I + (rho * R)
    @. dS[:, 1] += (-new_Svax * nu + new_Swane * psi)
    @. dS[:, 2] += (new_Svax * nu - new_Swane * psi)

    # change in exposed
    @. dE = new_I - (sigma * E)

    # calculate exposed to Is and Ia
    E_sigma = (sigma * E) * p_sigma
    E_inv_sigma = (sigma * E) * (1 - p_sigma)

    # change in infectious symptomatic
    @. dIs = E_sigma - ((gamma_Is .+ eta) .* Is)

    # change in infectious asymptomatic
    @. dIa = E_inv_sigma - (gamma_Ia * Ia)

    # change in hospitalised
    @. dH = (eta .* Is) - ((gamma_H + omega) .* H)

    # change in recovered
    @. dR = (gamma_Ia * Ia) + (gamma_Is * Is) +
            (gamma_H .* H) - (rho * R)

    @. dR[:, 1] += (-new_Rvax * nu + new_Rwane * psi)
    @. dR[:, 2] += (new_Rvax * nu - new_Rwane * psi)

    # change in dead
    @. dD = omega .* H
end

end


