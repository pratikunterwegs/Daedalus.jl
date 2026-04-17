```@meta
CurrentModule = Daedalus
```

# Vaccination Campaigns

Vaccination is modeled as a first-class **Event** type in Daedalus.jl.
`Vaccination` events generate ODE callbacks that modify epidemiological parameters during the simulation.

## Basic Vaccination Model

A `Vaccination` event tracks and controls:

- **Vaccination Rate (`rate`)**: Fraction of the population vaccinated per day
- **Uptake Limit (`uptake_limit`)**: Maximum fraction of the population willing to be vaccinated (0-1)
- **Efficacy (`efficacy`)**: Vaccine effectiveness against transmission (0-1)
- **Waning Period (`waning_period`)**: Days until vaccine-induced immunity decays
- **Campaign Timing (`start_time`)**: Day the vaccination campaign begins

## Simple Vaccination Campaign

Create a basic vaccination campaign and observe its effect on disease dynamics:

```@example basic_vaccination
using Daedalus
using Plots

# Define pathogen and vaccination campaign
infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")

vaccination = Daedalus.DaedalusStructs.Vaccination(
    100.0,    # Start at day 100
    0.001,          # Vaccinate 0.1% of population per day
    0.8,    # Up to 80% will accept vaccine
    0.9,        # 90% effective against transmission
    365.0  # Protection lasts ~1 year
)

# Run simulation with vaccination
data_with_vax = daedalus("Australia", infection, 
    vaccination=vaccination, time_end=500.0);

# Run baseline without vaccination for comparison
data_no_vax = daedalus("Australia", infection, time_end=500.0);

# Extract incidence and plot
times = Daedalus.Outputs.get_times(data_with_vax)
incidence_vax = Daedalus.Outputs.get_incidence(data_with_vax)
incidence_no_vax = Daedalus.Outputs.get_incidence(data_no_vax)

plot(times, incidence_vax, label="With vaccination", linewidth=2)
plot!(times, incidence_no_vax, label="No vaccination", linewidth=2)
xlabel!("Time (days)")
ylabel!("Daily incidence")
title!("Impact of Vaccination on Disease Incidence")
```

## Combining Vaccination with NPIs

Vaccination campaigns can be combined with non-pharmaceutical interventions for realistic policy scenarios. For example, vaccination might run during a period when transmission is being reduced through social distancing:

```@example vaccination_plus_npi
using Daedalus
using Plots

infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")

# Create a timed NPI: reduce transmission by 30% from day 50-200
npi_effect = Daedalus.DaedalusStructs.ParamEffect(
    :beta,
    x -> x .* 0.7,  # 30% reduction
    x -> x ./ 0.7,  # reset function
    Daedalus.DaedalusStructs.TimeTrigger(50.0),
    Daedalus.DaedalusStructs.TimeTrigger(200.0)
)
npi = Daedalus.DaedalusStructs.Npi([npi_effect])

# Vaccination starts after lockdown begins
vaccination = Daedalus.DaedalusStructs.Vaccination(
    start_time = 80.0,
    rate = 0.0015,      # 0.15% per day during lockdown
    uptake_limit = 0.75,
    efficacy = 0.85,
    waning_period = 365.0
)

# Run with both interventions
result = daedalus("Australia", infection, npi=npi, vaccination=vaccination, time_end=400.0)

times = Daedalus.Outputs.get_times(result)
incidence = Daedalus.Outputs.get_incidence(result)
rt = Daedalus.Outputs.get_values(result, "Rt", 1)

# Plot both epidemiological indicators
p1 = plot(times, incidence, label="Incidence", ylabel="Daily cases", legend=:topright)
p2 = plot(times, rt, label="Rt", ylabel="Reproduction number", legend=:topright)
plot(p1, p2, layout=(2,1), xlabel="Time (days)")
```

## Multiple Vaccination Phases

Model realistic multi-phase vaccination rollouts (e.g., healthcare workers first, then general population):

```@example multi_phase_vaccination
using Daedalus

infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")

# Phase 1: Slow healthcare worker vaccination (day 50-100)
vaccination_phase1 = Daedalus.DaedalusStructs.Vaccination(
    start_time = 50.0,
    rate = 0.0005,    # 0.05% per day (slow, targeted)
    uptake_limit = 0.2, # Only 20% of population (HCWs)
    efficacy = 0.95,
    waning_period = 365.0
)

# Phase 2: Accelerated general population vaccination (day 100+)
vaccination_phase2 = Daedalus.DaedalusStructs.Vaccination(
    start_time = 100.0,
    rate = 0.002,     # 0.2% per day (faster)
    uptake_limit = 0.7, # Up to 70% total
    efficacy = 0.85,  # Slightly lower efficacy in general population
    waning_period = 365.0
)

# Run with both vaccination phases
result = daedalus("Australia", infection, 
                 vaccination=[vaccination_phase1, vaccination_phase2],
                 time_end=500.0)

# Access the vaccination events from the result
vax_events = Daedalus.DaedalusOutput.get_vaccination(result)
println("Number of vaccination phases: $(length(vax_events))")
for (i, vax) in enumerate(vax_events)
    println("Phase $i: starts at day $(vax.start_time), efficacy $(vax.efficacy)")
end
```

## Understanding Vaccination Dynamics

During a vaccination campaign, the model:

1. **Tracks Coverage**: Accumulates vaccinated individuals each day until reaching the uptake limit
2. **Modifies Transmission**: Vaccinated populations experience reduced transmission (scaled by $1 - \text{efficacy}$)
3. **Applies Waning**: Immunity decreases over time according to the waning period

The vaccination rate is dynamic: it decreases automatically as coverage approaches the uptake limit, reflecting realistic saturation effects.

## Extracting Vaccination Metrics

After a simulation with vaccination, extract key metrics for analysis:

```julia
# Get all vaccination events used
vax_events = Daedalus.DaedalusOutput.get_vaccination(result)

# Check if a specific vaccination campaign was active
if !isempty(vax_events)
    vax = vax_events[1]
    @printf "Vaccination efficacy: %.1f%%\n" vax.efficacy * 100
    @printf "Uptake limit: %.1f%%\n" vax.uptake_limit * 100
    @printf "Campaign started: day %.0f\n" vax.start_time
end

# Get incidence and compare to baseline
incidence = Daedalus.Outputs.get_incidence(result)
total_infections = sum(incidence)
println("Total infections with vaccination: $(Int(round(total_infections)))")
```
