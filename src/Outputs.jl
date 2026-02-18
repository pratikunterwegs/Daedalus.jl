
module Outputs

using ..Constants

export get_values, get_times

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
output = daedalus(time_end=365.0)

# Daily hospitalizations
daily_H = get_values(output, "H", 1)

# Quarterly totals (90-day bins)
quarterly_H = get_values(output, "H", 90)

# Daily effective reproduction number
rt = get_values(output, "Rt", 1)
```
"""
function get_values(output, comp::String, timebin::Int = 90, strata::Union{Nothing,Vector{Int},UnitRange} = nothing)
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

end
