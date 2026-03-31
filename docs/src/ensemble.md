```@meta
CurrentModule = Daedalus
```

# Running multiple R0 values in parallel

The `daedalus` function supports running multiple reproduction number (R0) values in a single call using **vector dispatch**. This is useful for sensitivity analyses, parameter sweeps, or generating epidemic curves across a range of transmission scenarios.

**NOTE**: This section is work in progress. This feature uses _DifferentialEquations.jl_ [functionality for working with ensemble models](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/).
Better functionality for working with ensemble outputs, using the inbuilt DE.jl functionality, will be added later.

## Basic ensemble run

Pass a `Vector{InfectionData}` with different R0 values:

```@example ensemble_basic
using Daedalus
using Plots

# Run model with three different R0 values
covid_delta = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
infections = [
    covid_delta, deepcopy(covid_delta), deepcopy(covid_delta)
]
infections[1].r0 = 1.5
infections[2].r0 = 2.5
infections[3].r0 = 3.5
results = daedalus("Australia", infections, time_end=600.0);
```

## Extracting and comparing results

Each element in the results vector is a named tuple with the ODE solution and associated metadata:

```@example ensemble_basic
# Extract solutions for each R0 value
times = Daedalus.Outputs.get_times(results[1]);

# Compare exposed cases across R0 values
plot(legend=:topleft)
for (i, result) in enumerate(results)
    exposed = Daedalus.Outputs.get_values(result, "E", 1)
    plot!(times, exposed, label="R0 = $(result.r0)")
end
xlabel!("Time (days)")
ylabel!("Exposed individuals")
title!("Epidemic dynamics across R0 values")
```

## Comparing peak hospitalizations

Extract and compare peak hospitalization across different transmission scenarios:

```@example ensemble_basic
# get hospitalization for each R0 value and plot
# NOTE: needs better ensemble output handling functionality

# Display comparison
plot(legend=:topright)
for (i, result) in enumerate(results)
    hosp = Daedalus.Outputs.get_values(result, "H", 1)
    plot!(times, hosp, label="R0 = $(result.r0)")
end

xlabel!("Time (days)")
ylabel!("# Total hospital demand")
title!("Hospital demand by R0")
```

## Effective reproduction number (Rt) across scenarios

Compare how the effective reproduction number evolves under different baseline transmission rates:

```@example ensemble_basic
# Extract Rt for each R0 scenario
plot(legend=:topright)

for (i, result) in enumerate(results)
    rt_vals = Daedalus.Outputs.get_values(result, "Rt", 1)
    plot!(times, rt_vals, label="R0 = $(result.r0)")
end

hline!([1.0], linestyle=:dash, color=:black, label="Critical threshold")
xlabel!("Time (days)")
ylabel!("Effective reproduction number (Rt)")
title!("Rt trajectories for different R0 values")
```

## Using ensemble runs with NPIs

Non-pharmaceutical interventions (NPIs) are applied consistently across all R0 values in an ensemble run:

```@example ensemble_npi
using Daedalus

# Define a reactive NPI: reduce transmission by 50% when hospitalizations exceed 10,000
effect = Daedalus.DaedalusStructs.ParamEffect(
    :beta,
    x -> x .* 0.5;  # 50% reduction
    on = ("H", 10000.0),
    off = ("Rt", 1.0)
)
npi = Daedalus.DaedalusStructs.Npi([effect])

# Run ensemble with same NPI applied to all R0 values
r0_range = [1.5, 2.0, 2.5, 3.0, 3.5]
infections_npi = [Daedalus.DataLoader.get_pathogen("sars-cov-2 delta") for _ in r0_range]
for (i, r0) in enumerate(r0_range)
    infections_npi[i].r0 = r0
end
results_with_npi = daedalus("Australia", infections_npi, npi=npi, time_end=600.0);
```

## Large-scale parameter sweeps

For broader sensitivity analyses, you can define finer-grained R0 ranges:

```@example ensemble_basic
# Fine-grained sweep from 1.0 to 5.0 in 0.5 increments
r0_sweep = collect(1.0:0.5:5.0)
infections_sweep = [Daedalus.DataLoader.get_pathogen("sars-cov-2 delta") for _ in r0_sweep]
for (i, r0) in enumerate(r0_sweep)
    infections_sweep[i].r0 = r0
end
results_sweep = daedalus("Australia", infections_sweep, time_end=300.0);

# NOTE: add plot of infectious or exposed
```

## Accessing raw solution objects

The underlying ODE solution for each R0 value is accessible via the `sol` field, allowing direct access to state variables at any timepoint:

```@example ensemble_basic
# Get first result
result = results[1]

# Access solution at specific time points
state_at_day_100 = result.sol(100.0)
state_at_day_200 = result.sol(200.0)
```
