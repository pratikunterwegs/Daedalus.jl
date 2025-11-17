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

## Reactive events

```@example benchmarking_reactive_event
# benchmark for an RTM-exercise run of 100 days with two timed events
using Daedalus
using BenchmarkTools

npi = Daedalus.DaedalusStructs.Npi(20000.0, (coef=0.7,));

@benchmark daedalus(r0=5.0, npi=npi, time_end=600.0)
```

```@example benchmarking_reactive_event
# shorter duration
@benchmark daedalus(r0=5.0, npi=npi, time_end=100.0)
```
