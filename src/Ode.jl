
module Ode
export daedalus_ode!

"""
    daedalus_ode!(du, u, p, t)

The ODE system for the DAEDALUS model. This function is intended to be called
    internally from `daedalus`.
"""
function daedalus_ode!(du, u, p, t)
    i_age_groups = 1:4
    i_working_age = 3
    i_econ_groups = 5:49

    # du auto-magically takes the type of u (?)
    # each element of the tuple is one of the required params
    contacts, cw, beta, sigma, p_sigma, epsilon,
    rho, eta, omega, gamma_Ia, gamma_Is, gamma_H, nu, psi,
    t_vax, switch = p

    # view the values of each compartment per age group
    # rows represent age groups, epi compartments are columns
    S = @view u[:, 1]
    E = @view u[:, 2]
    Is = @view u[:, 3]
    Ia = @view u[:, 4]
    H = @view u[:, 5]
    R = @view u[:, 6]
    D = @view u[:, 7]
    V = @view u[:, 8]

    # calculate new infections and re-infections
    community_infectious = sum(Is .+ Ia * epsilon, dims=2)
    foi = beta * contacts * community_infectious

    # NOTE: element-wise multiplication
    new_I = S .* foi

    # workplace
    new_I_work = S[i_econ_groups] .* cw .*
                 (Is[i_econ_groups] .+ (Ia[i_econ_groups] * epsilon))

    new_I[i_econ_groups] += new_I_work

    # views to the change array slice
    dS = @view du[:, 1]
    dE = @view du[:, 2]
    dIs = @view du[:, 3]
    dIa = @view du[:, 4]
    dH = @view du[:, 5]
    dR = @view du[:, 6]
    dD = @view du[:, 7]
    dV = @view du[:, 8]

    # new vaccinations from susceptibles
    # new_Svax = @view S[:, 1]
    # new_Swane = @view S[:, 2]

    # new vaccinations from recovered
    # new_Rvax = @view R[:, 1]
    # new_Rwane = @view R[:, 2]

    # calculate change in compartment size and update the change matrix dU
    # note the use of @. for broadcasting, equivalent to .=
    # change in susceptibles
    @. dS = -new_I + (rho * R) - (nu * switch * S) + (psi * V)
    # @. dS[:, 1] += (-new_Svax * (nu * switch) + new_Swane * psi)
    # @. dS[:, 2] += (new_Svax * (nu * switch) - new_Swane * psi)

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
            (gamma_H .* H) - (rho + (nu * switch) * R)

    # @. dR[:, 1] += (-new_Rvax * (nu * switch) + new_Rwane * psi)
    # @. dR[:, 2] += (new_Rvax * (nu * switch) - new_Rwane * psi)

    # change in dead
    @. dD = omega .* H

    @. dV = nu * (S + R) - (psi * V)
end

end


