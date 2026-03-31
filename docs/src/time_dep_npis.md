
```@meta
CurrentModule = Daedalus
```

# Time-dependent interventions

`TimedNpi` interventions are responsive only to time, not epidemic state. Use these when intervention timing is fixed (e.g., legislative lockdowns or pre-planned public health campaigns).

### Single-Phase Example

```@example timed_npi
using Daedalus
using Plots

# Simple lockdown: reduce transmission by 70% from day 15 to day 45
timed_npi = Daedalus.DaedalusStructs.TimedNpi(
    30.0, # start time (days)
    90.0, # end time (days)
    0.3, # coefficient (0.3 = 70% reduction)
    "lockdown"
)

# Create infection data
infection_timed = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")

data = daedalus("Australia", "sars-cov-2 delta", npi=timed_npi, time_end=300.0);
data_default = daedalus("Australia", "sars-cov-2 delta", time_end=300.0);

# plot exposed
exposed = Daedalus.Outputs.get_values(data, "E", 1)
exposed_default = Daedalus.Outputs.get_values(data_default, "E", 1)
times = Daedalus.Outputs.get_times(data)

plot(times, exposed, label="Exposed w/ NPI")
xlabel!("Time (days)")
ylabel!("Exposed count")
```

### Multi-phased timed-NPI

Chain multiple intervention phases (e.g., lockdown → relaxation → second wave control):

```@example timed_npi_multi
using Daedalus
using Plots

# Define three intervention phases
timed_npi = Daedalus.DaedalusStructs.TimedNpi(
    [60.0, 180.0, 300.0],      # start times
    [120.0, 200.0, 365.0],     # end times
    [0.3, 0.7, 0.5],           # coefficients (70%, 30%, 50% reduction)
    "three_phase_lockdown"
)

data = daedalus("Australia", "sars-cov-2 delta", npi=timed_npi, time_end=400.0);

# Get exposed
exposed = Daedalus.Outputs.get_values(data, "E", 1)
times = Daedalus.Outputs.get_times(data)

plot(times, exposed, label="Exposed w/ NPI")
plot!(times, exposed_default, label="Exposed w/o NPI")
xlabel!("Time (days)")
ylabel!("Exposed count")
title!("Three-Phase Intervention Strategy")
```
