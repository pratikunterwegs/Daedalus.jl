```@meta
CurrentModule = Daedalus
```

This section shows how to run the Daedalus model.

```@example basic_daedalus
using Daedalus
using Plots

# pump up r0 to get peak within 50 days
data = daedalus(r0=5.0);

# plot exposed group
iExposed = Daedalus.Constants.get_indices("E")
iHosp = Daedalus.Constants.get_indices("H")
exposed = [sum(x[iExposed]) for x in data.sol.u]
hosp = [sum(x[iHosp]) for x in data.sol.u]

# plot the output to see lag in hospitalisations
plot(data.sol.t, exposed, label="exposed")
plot!(data.sol.t, hosp, label="hosp")
xlabel!("Time (days)")
ylabel!("# individuals")
```

```@example basic_daedalus
# plot recorded Rt
iRt = Daedalus.Constants.get_indices("Rt")
plot(data.sol.t, [x[iRt] for x in data.sol.u], label="Rt")
xlabel!("Time (days)")
ylabel!("Rt")
```

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
