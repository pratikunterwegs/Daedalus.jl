```@meta
CurrentModule = Daedalus
```

This module show how to use the `Npi` type to specify reactive events.

## Reactive NPIs

Updated to use the general `Npi` class rather than specialised `TimedNpi` and `ReactiveNpi` classes.

```@example reactive_npi
using Daedalus
using Plots

# create an Npi with the easy to use constructor
hosp_threshold=20000.0
a = Daedalus.DaedalusStructs.Npi(hosp_threshold, (coef=0.7,));

data = daedalus(r0=3.0, npi=a, time_end=600.0);
data_default = daedalus(r0=3.0, time_end=600.0);

# plot new hosp over time
iHosp = Daedalus.Constants.get_indices("H")
hosp = [sum(x[iHosp]) for x in data.sol(0:600)]
hosp_default = [sum(x[iHosp]) for x in data_default.sol(0:600)]

iRt = Daedalus.Constants.get_indices("Rt")
rt = [x[iRt] for x in data.sol(0:600)]
rt_default = [x[iRt] for x in data_default.sol(0:600)]
times = unique(data.sol.t)

plot(times, hosp, label="Hosp w/ NPI")
plot!(times, hosp_default, label="Hosp w/o NPI")
xlabel!("Time (days)")
ylabel!("Hosp occupancy")
```

Plot $R_t$ to see more differences.

```@example reactive_npi
plot(times, rt, label="Rt w/ NPI")
plot!(times, rt_default, label="Rt w/o NPI")
xlabel!("Time (days)")
ylabel!("Rt")
```

## Timed NPIs

Timed NPIs are responsive only to time and not state.
These were translated from the implementation in the R package _daedalus_ using Claude Sonnet 4.5 as a way of trying out Claude.

```@example timed_npi
using Daedalus
using Plots

# define a timed npi giving the start times, end times, and beta scaling
timed_npi = Daedalus.DaedalusStructs.TimedNpi(
    [60.0, 180.0, 300.0],
    [120.0, 200.0, 365.0],
    [0.3, 0.7, 0.5],
    "three_phase_lockdown"
)

data = daedalus(r0=3.0, npi=timed_npi, time_end=600.0);

# get exposed
iExpo = Daedalus.Constants.get_indices("E")
exposed = [sum(x[iExpo]) for x in data.sol(0:600)]
times = unique(data.sol.t)

plot(times, exposed, label="Incidence w/ NPI")
xlabel!("Time (days)")
ylabel!("Incidence")
```
