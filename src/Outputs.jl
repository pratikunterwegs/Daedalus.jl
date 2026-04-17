
module Outputs

using ..Constants
using DataFrames

export get_values, get_times, get_incidence, get_epidemic_summary, get_life_years_lost

"""
    get_values(output, comp::String, timebin::Int=90)

Extract values for a given compartment from a `daedalus()` model output, with
optional aggregation into time bins.

Assumes the model was solved with daily savepoints (integer timesteps from 0).

# Arguments
- `output`: Named tuple returned by `daedalus()`, with field `sol` (the ODE solution).
- `comp`: Compartment name. One of: `"S"`, `"E"`, `"Is"`, `"Ia"`, `"H"`, `"R"`,
  `"D"`, or `"Rt"`.
- `timebin`: Number of days per aggregation bin. Use `1` to return unaggregated
  daily values (default: `90`).

# Returns

A `Vector{Float64}` where each element is the sum of the compartment values across
all age and vaccination groups:

- `timebin == 1`: one value per day, from day 0 to day `tmax` (length `tmax + 1`).
- `timebin > 1`: one value per non-overlapping window of `timebin` days. The final
  bin covers the remaining days and may be shorter than `timebin`.

# Examples

```julia
output = daedalus(country="Australia", time_end=365.0)

# Daily hospitalizations
daily_H = get_values(output, "H", 1)

# Quarterly totals (90-day bins)
quarterly_H = get_values(output, "H", 90)

# Daily effective reproduction number
rt = get_values(output, "Rt", 1)
```
"""
function get_values(output, comp::String, timebin::Int = 90,
        strata::Union{Nothing, Vector{Int}, UnitRange} = nothing)
    compidx = get_indices(comp)
    tmax = maximum(output.sol.t)
    times = output.sol.t
    u = output.sol.u

    if isnothing(strata)
        compidx = get_indices(comp)
    else
        compidx = get_indices(comp, strata)
    end

    # Extract the summed compartment value at each integer timestep
    n_steps = round(Int, tmax) + 1
    vals = Vector{Float64}(undef, n_steps)
    for (i, t) in enumerate(0.0:1.0:tmax)
        tidx = findfirst(x -> isapprox(x, t), times)
        vals[i] = sum(u[tidx][compidx])
    end

    timebin == 1 && return vals

    # Aggregate into non-overlapping bins using cumulative sums.
    # bin_ends[b] is the exclusive upper boundary of bin b (1-indexed into vals),
    # so bin b covers vals[(b-1)*timebin+1 : bin_ends[b]].
    tmax_int = round(Int, tmax)
    bin_ends = unique([collect(timebin:timebin:tmax_int); tmax_int]) .+ 1
    cumulative = cumsum(vals)
    binned = cumulative[bin_ends]

    return [binned[1]; diff(binned)]
end

"""
    get_times(output)

Return the unique timepoints at which the model state was saved.

# Arguments
- `output`: Named tuple returned by `daedalus()`, with field `sol` (the ODE solution).

# Returns
A sorted `Vector{Float64}` of unique save times (e.g. `[0.0, 1.0, 2.0, ..., tmax]`
for a daily-increment run).
"""
function get_times(output)
    return unique(output.sol.t)
end

"""
    get_incidence(output::DaedalusOutput)

Compute daily incidence (new infections, hospitalisations, deaths) from a simulation.

Returns a DataFrame with columns: time, measure, value
"""
function get_incidence(output)
    sol = output.sol
    tmax = maximum(sol.t)
    times = 0:1:Int(round(tmax))

    # Extract daily values
    infections = Float64[]
    hosp = Float64[]
    deaths = Float64[]

    for t in times
        tidx = findfirst(x -> isapprox(x, t), sol.t)
        if !isnothing(tidx)
            u = sol.u[tidx]
            inf_val = sum(u[get_indices("new_inf")])
            push!(infections, inf_val)
            h_val = sum(u[get_indices("new_hosp")])
            push!(hosp, h_val)
            d_val = sum(u[get_indices("D")])
            push!(deaths, d_val)
        end
    end

    # Compute daily differences
    inf_daily = diff(infections)
    hosp_daily = diff(hosp)
    death_daily = diff(deaths)

    # Build DataFrame
    n_times = length(times) - 1
    times_vec = Vector(times[2:end])
    measures = vcat([repeat([m], n_times)
                     for m in ["infections", "hospitalisations", "deaths"]]...)
    values = vcat(inf_daily, hosp_daily, death_daily)
    times_vec_rep = vcat([times_vec for _ in 1:3]...)

    return DataFrame(time = times_vec_rep, measure = measures, value = values)
end

"""
    get_epidemic_summary(output::DaedalusOutput)

Compute total epidemic summary (infections, hospitalisations, deaths).

Returns a DataFrame with columns: measure, value
"""
function get_epidemic_summary(output)
    sol = output.sol
    u_final = sol.u[end]

    # Final cumulative values
    total_inf = sum(u_final[get_indices("E")]) + sum(u_final[get_indices("Is")]) +
                sum(u_final[get_indices("Ia")]) + sum(u_final[get_indices("R")]) +
                sum(u_final[get_indices("D")])
    total_hosp = sum(u_final[get_indices("H")]) + sum(u_final[get_indices("R")]) +
                 sum(u_final[get_indices("D")])
    total_deaths = sum(u_final[get_indices("D")])

    return DataFrame(
        measure = ["infections", "hospitalisations", "deaths"],
        value = [total_inf, total_hosp, total_deaths]
    )
end

"""
    get_life_years_lost(output::DaedalusOutput)

Compute life-years lost due to deaths.

Returns a NamedTuple with fields: total, by_age
"""
function get_life_years_lost(output)
    sol = output.sol
    u_final = sol.u[end]
    u_initial = sol.u[1]

    # Deaths per 4 age groups (ages 0-4, 5-19, 20-64, 65+)
    d_idx = get_indices("D")
    deaths_by_age = [sum(u_final[d_idx[get_age_group_indices(ag)]]) -
                     sum(u_initial[d_idx[get_age_group_indices(ag)]])
                     for ag in 1:Constants.N_AGE_GROUPS]

    lyl_by_age = deaths_by_age .* output.country.life_expectancy
    lyl_total = sum(lyl_by_age)

    return (total = lyl_total, by_age = lyl_by_age)
end

end
