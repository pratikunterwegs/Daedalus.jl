```@meta
CurrentModule = Daedalus
```

# Implementing interventions using Daedalus Event types

Daedalus.jl uses a unified **Event** framework for managing interventions and campaigns during epidemics.
Both non-pharmaceutical interventions (NPIs) and vaccination campaigns are modeled as Event types, making them composable and extensible.

## Event Types

All events inherit from an abstract `Event` base type. Currently, two concrete event types are implemented:

### Npi (Non-Pharmaceutical Interventions)

`Npi` represents policy-based or behavioral interventions that modify epidemiological parameters. Each `Npi` contains one or more `ParamEffect` objects, each specifying:
- Which parameter to modify (e.g., `:beta` for transmission)
- How to modify it (transformation function)
- When to activate and deactivate (via `Trigger` conditions)

### Vaccination

`Vaccination` represents a vaccination campaign that dynamically modifies vaccination rates (`:nu`) and immunity waning rates (`:psi`) throughout the simulation.

## Event Composition

Events can be freely combined in a single simulation:

```@example combined_events
using Daedalus
using Plots

# Create a time-limited NPI: reduce transmission by 40% from day 30-90
infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
effect = Daedalus.DaedalusStructs.ParamEffect(
    :beta, x -> x .* 0.6, x -> x ./ 0.6,
    Daedalus.DaedalusStructs.TimeTrigger(30.0),
    Daedalus.DaedalusStructs.TimeTrigger(90.0)
)
npi = Daedalus.DaedalusStructs.Npi([effect])

# Create a vaccination campaign
vaccination = Daedalus.DaedalusStructs.Vaccination(
    50.0,    # begins at day 50
    0.001,         # 0.1% of population per day
    0.7,   # max 70% willing to be vaccinated
    0.9,       # 90% effective
    365.0 # immunity wanes after 1 year
)

# Run with both NPI and vaccination
result = daedalus("Australia", infection, npi=npi, 
    vaccination=vaccination, time_end=300.0);

times = Daedalus.Outputs.get_times(result);
exposed = Daedalus.Outputs.get_values(result, "E", 1);

plot(times, exposed, label="Exposed", linewidth=2)
xlabel!("Time (days)")
ylabel!("# Exposed")
title!("Epidemic with NPI and Vaccination")
```

## Extracting events from results

The `DaedalusOutput` object contains all events that were active during the simulation in its `.events` field. Use accessor functions to retrieve specific event types:

```julia
# Get all NPIs from the result
npis = Daedalus.DaedalusOutput.get_npi(result);

# Get all vaccination campaigns
vaccines = Daedalus.DaedalusOutput.get_vaccination(result)
```
