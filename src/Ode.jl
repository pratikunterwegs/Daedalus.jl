
module Ode
export daedalus_ode!

using ..Constants
using ..DaedalusStructs
using ..Helpers

using LinearAlgebra
using StaticArrays

"""
    daedalus_ode!(du, u, p, t)

Compute derivatives for the DAEDALUS epidemic ODE system.

The core compartmental model with 9 compartments (7 primary: S, E, Is, Ia, H, R, D;
2 data: newI, newH) across 49 demographic groups (4 age + 45 economic sectors)
and 2 vaccine strata. Implements force-of-infection, disease progression, vaccination,
and death.

# Arguments
- `du::Array`: Pre-allocated output array (derivatives), same shape as `u`
- `u::Array`: State vector (length = N_TOTAL_GROUPS × N_COMPARTMENTS × N_VACCINE_STRATA + new_vax + Rt)
- `p::Params`: Parameter struct containing contact matrices, rates, and NGM
- `t::Number`: Current simulation time

# State layout
The state array is reshaped internally as (49, 9, 2) where:
- Dimension 1: Demographic groups (4 age + 45 economic sectors)
- Dimension 2: Compartments (S, E, Is, Ia, H, R, D, newI, newH) — last 2 are cumulative data
- Dimension 3: Vaccine strata (unvaccinated, vaccinated)
- new_vax block: Cumulative new vaccinations per group
- Final element: Effective reproduction number Rt (updated by callbacks)

# Details
This function is called internally by OrdinaryDiffEq's ODE solvers. It computes
force-of-infection from the contact matrix and current infectious compartments,
then updates all compartments according to the SEIR-like disease dynamics.
"""
function daedalus_ode!(du::Array, u::Array, p::Params, t::Number)
    # du auto-magically takes the dims and type of u (?)
    # each element of the tuple is one of the required params
    # time dependent vaccination
    nu_eff = p.nu

    # view the values of each compartment per age group
    # rows represent age groups, epi compartments are columns
    U = @view u[1:p.size_state]
    U = reshape(U, (N_TOTAL_GROUPS, N_COMPARTMENTS, N_VACCINE_STRATA))

    dU = @view du[1:p.size_state]
    dU = reshape(dU, (N_TOTAL_GROUPS, N_COMPARTMENTS, N_VACCINE_STRATA))

    # using views does not seem to affect performance greatly
    S = @view U[:, iS, :]
    E = @view U[:, iE, :]
    Is = @view U[:, iIs, :]
    Ia = @view U[:, iIa, :]
    H = @view U[:, iH, :]
    R = @view U[:, iR, :]

    # calculate new infections and re-infections
    community_infectious = sum(Is .+ Ia .* p.epsilon, dims = 2)

    # NOTE: this can probably be cleaned up
    # NOTE: this is inefficient and should be refactored so that p.cm_temp
    # is modified only when NPIs fire. this is a proof-of-concept.
    cm_scaling = ones(p.settings) # NOTE: placeholder for future operations
    weighted_slice_sum!(p.contacts, cm_scaling, p.cm_temp)

    foi = p.beta * p.cm_temp * community_infectious

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
    @. dH = (p.eta .* Is) - ((p.gamma_H + p.omega) .* H)

    # change in recovered
    @. dR = (p.gamma_Ia * Ia) + (p.gamma_Is * Is) +
            (p.gamma_H .* H) - (p.rho * R)

    @. dR[:, 1] += (-new_Rvax * nu_eff + new_Rwane * p.psi)
    @. dR[:, 2] += (new_Rvax * nu_eff - new_Rwane * p.psi)

    # change in dead
    @. dD = p.omega .* H

    # Track cumulative new infections and new hospitalisations (data compartments)
    dnewI = @view dU[:, inewI, :]
    @. dnewI = new_I
    dnewH = @view dU[:, inewH, :]
    @. dnewH = p.eta .* Is

    # Track cumulative new vaccinations per group (matches C++ t_new_vax)
    new_vax_view = @view du[(p.size_state + 1):(p.size_state + N_TOTAL_GROUPS)]
    @. new_vax_view = nu_eff * new_Svax + nu_eff * new_Rvax

    # Rt is updated by callbacks; keep its derivative zero between updates
    du[end] = 0.0
end

end
