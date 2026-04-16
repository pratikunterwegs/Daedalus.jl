module Costs

using ..Constants

export get_costs

"""
    get_costs(output::DaedalusOutput)

Compute economic and life-related costs from a completed simulation.

Returns a NamedTuple with fields:
- total_cost: sum of economic and life value costs
- economic_costs: NamedTuple with economic loss breakdown
- life_value_lost: NamedTuple with (total, by_age)
- life_years_lost: NamedTuple with (total, by_age)
"""
function get_costs(output)
    sol = output.sol
    cd = output.country

    u_final = sol.u[end]
    u_initial = sol.u[1]

    # Deaths by 4 age groups (ages 0-4, 5-19, 20-64, 65+)
    # State layout: N_VACCINE_STRATA * N_COMPARTMENTS * N_TOTAL_GROUPS
    # D compartment is at index 7 in epi compartments
    d_idx = get_indices("D")

    deaths_by_age = [
        sum(u_final[d_idx[get_age_group_indices(ag)]]) - sum(u_initial[d_idx[get_age_group_indices(ag)]])
        for ag in 1:Constants.N_AGE_GROUPS
    ]

    # Life years lost = deaths * remaining life expectancy
    lyl_by_age = deaths_by_age .* cd.life_expectancy
    lyl_total = sum(lyl_by_age)

    # Life value lost = deaths * VSL by age
    lvl_by_age = deaths_by_age .* cd.vsl
    lvl_total = sum(lvl_by_age)

    # Simple economic cost approximation: worker absences due to infection
    # Is_sector * (gva/workers) over time
    is_idx = get_indices("Is")
    econ_loss = 0.0

    for t_step in eachindex(sol.t)
        u_t = sol.u[t_step]
        # Sum Is across working-age population (age groups 3, economic groups 1-45)
        for eg in 1:Constants.N_ECON_GROUPS
            state_idx = Constants.N_TOTAL_GROUPS * (Constants.iIs - 1) + Constants.N_AGE_GROUPS + eg
            if 1 <= state_idx <= length(u_t)
                is_val = u_t[state_idx]
                gva_per_worker = cd.gva[eg] / max(cd.workers[eg], 1)
                econ_loss += is_val * gva_per_worker
            end
        end
    end

    return (
        total_cost = econ_loss + lvl_total,
        economic_costs = (total = econ_loss,),
        life_value_lost = (total = lvl_total, by_age = lvl_by_age),
        life_years_lost = (total = lyl_total, by_age = lyl_by_age)
    )
end

end
