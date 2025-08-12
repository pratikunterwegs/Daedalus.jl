```@meta
CurrentModule = Daedalus
```

This section shows how to run the Daedalus model.

```@example basic_daedalus
using Daedalus
using Plots

# all arguments have appropriate defaults
data = daedalus()

# plot exposed group
plot(data, vars=(0, 50:99))
```

```@example basic_daedalus
# plot exposed and vaccinated group
# check that vaccinations begin at t = t_vax = 200.0
plot(data, vars=(0, 343:392))
```

## Intervention

```@example timed_npi
using Daedalus
using Plots

# make an NPI
params = (coef = 0.2,)
time_on = [30.0, 50.0]
time_off = [45.0, 65.0]

npi = Daedalus.DaedalusStructs.TimedNpi(
    params=params, time_on=time_on, time_off=time_off
)

data = daedalus(beta=0.1, npi=npi)

# plot output
plot(data, vars=(0, 50:99))
```

## State-dependent event

This example shows a reactive event triggered by a state: mortality rate for all groups increases by 1.6x when hospital capacity is breached.

```@example reactive_event
using Daedalus
using Plots

# make reactive NPI
idx = Daedalus.Constants.get_indices("H")
hosp_capacity = 10000.0

npi = Daedalus.DaedalusStructs.ReactiveNpi(
    params=(coef=1.6,),
    id_state_on=idx, id_state_off=idx,
    value_state_on=hosp_capacity, value_state_off=hosp_capacity-100.0
)

Daedalus.Events.get_coef(npi)

data = daedalus(beta = 0.05, npi=npi, time_end=200.0)
data_default = daedalus(beta = 0.05, time_end=200.0)

idx_deaths = Daedalus.Constants.get_indices("D")
idx_hosp = Daedalus.Constants.get_indices("H")

# plot new hosp over time
new_hosp = map((x) -> sum(x[idx_hosp]), data.u)
new_hosp_default = map((x) -> sum(x[idx_hosp]), data_default.u)

plot(new_hosp_default)
plot!(new_hosp, linecolor="orange")
hline!([hosp_capacity], linestyle=:dash)
```

```@example reactive_event
# plot new deaths over time
new_deaths = diff(map((x) -> sum(x[idx_deaths]), data.u))
new_deaths_default = diff(map((x) -> sum(x[idx_deaths]), data_default.u))

plot(new_deaths_default)
plot!(new_deaths, linecolor="orange")
```

## IPR-dependent NPI

```@example ipr_event
using Daedalus
using Plots

# make reactive NPI
idx_on = Daedalus.Constants.get_indices("H")
hosp_capacity = 10000.0

idx_off = 687 # manual
ipr_limit = 0.476 # manual value of gamma_Ia

npi = Daedalus.DaedalusStructs.ReactiveNpi(
    params=(coef=0.2,),
    id_state_on=idx_on, id_state_off=[idx_off],
    value_state_on=hosp_capacity, value_state_off=ipr_limit
)

Daedalus.Events.get_coef(npi)

data = daedalus(beta = 1.3 / 7.0, time_end=30.0, print=false)
data_default = daedalus(beta = 0.05, time_end=200.0)

idx_exposed = Daedalus.Constants.get_indices("E")

# plot new hosp over time
new_exposed = map((x) -> sum(x[idx_exposed]), data.u)
new_exposed_default = map((x) -> sum(x[idx_exposed]), data_default.u)

ipr = map((x) -> sum(x[(idx_off +1)]), data.u)
plot(ipr)

plot(new_exposed_default)
plot(new_exposed, linecolor="orange")
```

```@example ipr_event
# plot new deaths over time
ipr = diff(map((x) -> sum(x[idx_off]), data.u))

plot(ipr)
hline!([ipr_limit], linestyle=:dash)
```

## Rt example

```@example ipr_event
using Daedalus
using Plots

idx_off = 688 # manual

data = daedalus(beta = 1.3 / 7.0, time_end=30.0, print=false);

plot(data.saves.t, data.saves.saveval)

rt = map((x) -> sum(x[(idx_off)]), data.sol.u)
plot!([0.0; data.saves.t], rt, color="red")

# plot new hosp over time
new_exposed = map((x) -> sum(x[idx_exposed]), data.u)
plot(new_exposed, linecolor="orange")
vline!([21.38])
```