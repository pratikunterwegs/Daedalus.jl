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

## Effect of logging $R_t$

Logging $R_t$ in each timestep slows the simulation down by _a lot_.

```@example benchmarking_rt_logging
using Daedalus
using BenchmarkTools

# turn off Rt logging and compare with benchmark above
@benchmark daedalus(log_rt=false, time_end=600.0)
```

## Reactive events

Reactive events currently require logging $R_t$ so are not benchmarked with logging turned off.

```@example benchmarking_reactive_event
# benchmark for an RTM-exercise run of 100 days
using Daedalus
using BenchmarkTools

npi = Daedalus.DaedalusStructs.Npi(20000.0, (coef=0.7,));

@benchmark daedalus(r0=5.0, npi=npi, time_end=600.0)
```

```@example benchmarking_reactive_event
# shorter duration
@benchmark daedalus(r0=5.0, npi=npi, time_end=100.0)
```

## Timed events

```@example benchmarking_timed_event
# benchmark for a full run of 600 days
using Daedalus
using BenchmarkTools

timed_npi = Daedalus.DaedalusStructs.TimedNpi(
    [60.0, 180.0, 300.0],
    [120.0, 200.0, 365.0],
    [0.3, 0.7, 0.5],
    "three_phase_lockdown"
)

@benchmark daedalus(r0=3.0, npi=timed_npi, time_end=600.0)
```

Timed events do not need $R_t$ logging and can be benchmarked without it.

```@example benchmarking_timed_event
@benchmark daedalus(r0=3.0, npi=timed_npi, time_end=600.0, log_rt=false)
```
