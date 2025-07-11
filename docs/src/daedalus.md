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

```@example npi
using Daedalus
using Plots

# make an NPI
params = (coef = 0.2,)
time_on = [30.0, 50.0]
time_off = [45.0, 65.0]

npi = Daedalus.DaedalusStructs.Npi(
    params=params, time_on=time_on, time_off=time_off
)

data = daedalus(beta=0.1, npi=npi)

# plot output
plot(data, vars=(0, 50:99))
```

## Benchmarking

```@example benchmarking
# benchmark for a typical daedalus run of 600 days
using Daedalus
using BenchmarkTools

@benchmark daedalus()
```
