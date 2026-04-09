
```@meta
CurrentModule = Daedalus
```

# Time-dependent interventions

`TimedNpi` interventions are responsive only to time, not epidemic state.

### Single-Phase Example

```@example timed_npi
using Daedalus
using Plots

# Simple lockdown: reduce transmission by 70% from day 15 to day 45
timed_npi = Npi(
    [TimedEffect(
        :beta, x -> x .* 0.3, x -> x ./ 0.3,
        30.0, 90.0
    )]
)

data = daedalus("Australia", "sars-cov-2 delta", npi=timed_npi, time_end=300.0);
data_default = daedalus("Australia", "sars-cov-2 delta", time_end=300.0);

# plot exposed
exposed = Daedalus.Outputs.get_values(data, "E", 1)
exposed_default = Daedalus.Outputs.get_values(data_default, "E", 1)
times = Daedalus.Outputs.get_times(data)

plot(times, exposed, label="Exposed w/ NPI")
plot!(times, exposed_default, label = "Exposed w/o NPI")
xlabel!("Time (days)")
ylabel!("# exposed")
```

### Multi-phased timed-NPI

Chain multiple intervention phases (e.g., lockdown → relaxation → second wave control):

```@example timed_npi_multi
using Daedalus
using Plots

# Define three intervention phases
timed_npi = Npi(
    [TimedEffect(
        :beta, x -> x .* 0.8, x -> x ./ 0.8, 60.0, 120.0
    ),
    TimedEffect(
        :beta, x -> x .* 0.5, x -> x ./ 0.7, 130.0, 200.0
    ),
    TimedEffect(
        :beta, x -> x .* 0.8, x -> x ./ 0.5, 210.0, 365.0
    )]
)

data = daedalus("Australia", "sars-cov-2 delta", npi=timed_npi, time_end=400.0);
data_default = daedalus("Australia", "sars-cov-2 delta", time_end=400.0);

# Get exposed
exposed = Daedalus.Outputs.get_values(data, "E", 1)
exposed_default = Daedalus.Outputs.get_values(data_default, "E", 1)
times = Daedalus.Outputs.get_times(data)

plot(times, exposed, label="Exposed w/ NPI")
plot!(times, exposed_default, label="Exposed w/o NPI")
xlabel!("Time (days)")
ylabel!("# exposed")
title!("Three-Phase Intervention Strategy")
```
