```@meta
CurrentModule = Daedalus
```

# NPI Helper Functions

Creating complex interventions by directly constructing `Npi`, `ParamEffect`, and trigger objects can be verbose. The **NPI helpers** module provides convenience functions that reduce boilerplate while maintaining full expressiveness.

## Overview

Three helper functions cover the most common intervention patterns:

- `create_time_intervention()`: Single-parameter time-limited interventions
- `create_reactive_intervention()`: State-dependent interventions (threshold-triggered)
- `create_multi_intervention()`: Multi-parameter coordinated interventions

All return ready-to-use `Npi` objects compatible with `daedalus()`.

## Time-Limited Interventions

Use `create_time_intervention()` for policies active only during a specific time window (e.g., lockdowns with known start/end dates).

### Simple Transmission Reduction

```@example time_npi_basic
using Daedalus
using Plots

infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")

# Reduce transmission by 60% from day 30 to day 90
npi = Daedalus.NpiHelpers.create_time_intervention(
    :beta,                          # target parameter
    x -> x .* 0.4,                  # reduction: multiply by 0.4 (60% reduction)
    start = 30.0,
    end_time = 90.0
)

result = daedalus("Australia", infection, npi=npi, time_end=300.0)
result_baseline = daedalus("Australia", infection, time_end=300.0)

# Compare epidemic curves
times = Daedalus.Outputs.get_times(result)
incidence = Daedalus.Outputs.get_incidence(result)
incidence_baseline = Daedalus.Outputs.get_incidence(result_baseline)

plot(times, incidence, label="With lockdown (days 30-90)", linewidth=2)
plot!(times, incidence_baseline, label="Baseline", linewidth=2)
xlabel!("Time (days)")
ylabel!("Daily incidence")
title!("Impact of 60-day Lockdown")
legend!(:topright)
```

### Other Parameters

The same helper works for any parameter:

```julia
# Slow hospitalization recovery (simulate hospital stress)
npi = create_time_intervention(
    :gamma_H,
    x -> x .* 0.8,  # 20% slower recovery
    start = 100.0,
    end_time = 200.0
)

# Reduce contact rates (simulate behavior change)
npi = create_time_intervention(
    :omega,
    x -> x .* 0.6,  # 40% fewer contacts
    start = 50.0,
    end_time = 150.0
)
```

## State-Dependent Interventions

Use `create_reactive_intervention()` for policies triggered by epidemic thresholds (e.g., "activate lockdown if hospitalizations exceed 5000").

### Hospital Capacity Trigger

```@example reactive_npi_helper
using Daedalus
using Plots

infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")

# Reduce transmission when hospital capacity is strained
npi = Daedalus.NpiHelpers.create_reactive_intervention(
    :beta,
    x -> x .* 0.5,                 # 50% transmission reduction
    on_threshold = 8000.0,         # activate if H > 8000
    off_threshold = 4000.0,        # deactivate if H < 4000
    compartment = "H"              # monitor hospitalizations
)

result = daedalus("Australia", infection, npi=npi, time_end=300.0)

times = Daedalus.Outputs.get_times(result)
hospitalized = Daedalus.Outputs.get_values(result, "H", 1)
rt = Daedalus.Outputs.get_values(result, "Rt", 1)

# Show when NPI was active (when H was elevated)
p1 = plot(times, hospitalized, label="Hospitalizations", ylabel="# hospitalized",
         fill=(0, 0.3, :lightblue), legend=:topleft)
hline!([8000, 4000], label=["NPI on" "NPI off"], linestyle=:dash, color=:red)

p2 = plot(times, rt, label="Rt", ylabel="Reproduction number", legend=:topright)
plot(p1, p2, layout=(2,1), xlabel="Time (days)")
```

### Monitor Different Compartments

```julia
# Activate if exposed cases exceed a threshold
npi = create_reactive_intervention(
    :beta,
    x -> x .* 0.6,
    on_threshold = 50000.0,
    off_threshold = 20000.0,
    compartment = "E"  # monitor exposed, not hospitalized
)

# Activate if reproduction number gets too high
npi = create_reactive_intervention(
    :beta,
    x -> x .* 0.7,
    on_threshold = 1.5,
    off_threshold = 1.0,
    compartment = "Rt"  # monitor Rt directly
)
```

## Multi-Parameter Interventions

Use `create_multi_intervention()` for coordinated policies affecting multiple parameters simultaneously (e.g., school and workplace closures affecting both transmission and contact rates).

### Comprehensive Lockdown

```@example multi_param_helper
using Daedalus
using Plots

infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")

# Coordinated policy: reduce transmission and increase recovery awareness
reductions = Dict(
    :beta => x -> x .* 0.3,        # 70% transmission reduction
    :omega => x -> x .* 1.1        # 10% faster recovery (better healthcare)
)

npi = Daedalus.NpiHelpers.create_multi_intervention(
    reductions,
    start = 50.0,
    end_time = 150.0
)

result = daedalus("Australia", infection, npi=npi, time_end=300.0)

times = Daedalus.Outputs.get_times(result)
exposed = Daedalus.Outputs.get_values(result, "E", 1)
recovered = Daedalus.Outputs.get_values(result, "R", 1)

plot(times, exposed, label="Exposed", linewidth=2)
plot!(times, recovered, label="Recovered", linewidth=2)
xlabel!("Time (days)")
ylabel!("# individuals")
title!("Multi-Parameter Intervention (Lockdown + Healthcare Surge)")
```

### Examples with Different Parameters

```julia
# School closure: reduce student transmission and contact patterns
npi = create_multi_intervention(
    Dict(
        :beta => x -> x .* 0.6,    # reduced transmission in schools
        :omega => x -> x .* 0.8    # reduced overall contacts
    ),
    start = 30.0,
    end_time = 100.0
)

# Healthcare surge: boost hospital capacity and recovery
npi = create_multi_intervention(
    Dict(
        :gamma_H => x -> x .* 1.3,    # faster hospital discharge
        :eta => x -> x .* 1.1         # better treatment outcomes
    ),
    start = 80.0,
    end_time = 180.0
)
```

## Combining Multiple Interventions

Create complex intervention scenarios by building multiple NPIs separately and passing them as a vector:

```@example multi_npi_scenario
using Daedalus

infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")

# Phase 1: Initial lockdown (days 30-60)
npi1 = Daedalus.NpiHelpers.create_time_intervention(
    :beta, x -> x .* 0.3, start = 30.0, end_time = 60.0
)

# Phase 2: Partial reopening with reactive trigger (days 60+)
npi2 = Daedalus.NpiHelpers.create_reactive_intervention(
    :beta, x -> x .* 0.7,
    on_threshold = 5000.0, off_threshold = 2000.0,
    compartment = "H"
)

# Combine into a vector to pass to daedalus
npis = [npi1, npi2]

# Create a unified NPI by flattening all effects
combined_effects = vcat([e.effects for e in npis]...)
combined_npi = Daedalus.DaedalusStructs.Npi(combined_effects)

result = daedalus("Australia", infection, npi=combined_npi, time_end=300.0)
```

## When to Use Helpers vs. Direct Construction

| Use Case | Recommended Approach |
|----------|---------------------|
| Simple single-parameter, timed intervention | `create_time_intervention()` |
| Simple state-dependent trigger | `create_reactive_intervention()` |
| Coordinated multi-parameter response | `create_multi_intervention()` |
| Custom trigger logic or non-standard patterns | Direct `ParamEffect` + `Trigger` construction |
| Complex event sequences | Combine multiple helpers or build `Npi` manually |

The helper functions are designed to be intuitive for common scenarios while remaining flexible enough for most modeling needs.
