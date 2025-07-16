```@meta
CurrentModule = Daedalus
```

# Benchmarking

This section shows various benchmarks.

## Runs of different durations

```@example benchmarking_01
# benchmark for a typical daedalus run of 600 days
using Daedalus
using BenchmarkTools

@benchmark daedalus(time_end=600.0)
```

```@example benchmarking_02
# benchmark for an RTM-exercise run of 100 days
using Daedalus
using BenchmarkTools

@benchmark daedalus(time_end=100.0)
```

## Timed events

```@example benchmarking_timed_events
# benchmark for an RTM-exercise run of 100 days with a single timed event
using Daedalus
using BenchmarkTools

params = (coef = 0.2,)
time_on = [30.0]
time_off = [65.0]

npi = Daedalus.DaedalusStructs.TimedNpi(
    params=params, time_on=time_on, time_off=time_off
)

@benchmark daedalus(npi=npi, time_end=100.0)
```

```@example benchmarking_multi_timed_events
# benchmark for an RTM-exercise run of 100 days with two timed events
using Daedalus
using BenchmarkTools

params = (coef = 0.2,)
time_on = [30.0, 50.0]
time_off = [45.0, 65.0]

npi = Daedalus.DaedalusStructs.TimedNpi(
    params=params, time_on=time_on, time_off=time_off
)

@benchmark daedalus(npi=npi, time_end=100.0)
```

## State-dependent events

```@example benchmarking_reactive_event
# benchmark for an RTM-exercise run of 100 days with two timed events
using Daedalus
using BenchmarkTools

idx = Daedalus.Constants.get_indices("H")
hosp_capacity = 10000.0

npi = Daedalus.DaedalusStructs.ReactiveNpi(
    params=(coef=1.6,),
    id_state_on=idx, id_state_off=idx,
    value_state_on=hosp_capacity, value_state_off=hosp_capacity-100.0
)

@benchmark daedalus(npi=npi, time_end=100.0)
```
