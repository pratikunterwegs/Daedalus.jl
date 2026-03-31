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

infection_01 = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
infection_01.r0 = 2.5
@benchmark daedalus("Australia", infection_01, time_end=600.0)
```

```@example benchmarking_02
# benchmark for an RTM-exercise run of 100 days
using Daedalus
using BenchmarkTools

infection_02 = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
infection_02.r0 = 2.5
@benchmark daedalus("Australia", infection_02, time_end=100.0)
```

## Effect of logging $R_t$

Logging $R_t$ in each timestep slows the simulation down by _a lot_.

```@example benchmarking_rt_logging
using Daedalus
using BenchmarkTools

# turn off Rt logging and compare with benchmark above
infection_rt = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
infection_rt.r0 = 2.5
@benchmark daedalus("Australia", infection_rt, log_rt=false, time_end=600.0)
```

## Reactive events

Reactive events currently require logging $R_t$ so are not benchmarked with logging turned off.

```@example benchmarking_reactive_event
# benchmark for an RTM-exercise run with reactive NPI
using Daedalus
using BenchmarkTools

# Create a reactive NPI: reduce transmission when hospitalizations exceed 20,000
effect = Daedalus.DaedalusStructs.ParamEffect(
    :beta,
    x -> x .* 0.3;  # 70% reduction
    on = ("H", 20000.0),
    off = ("Rt", 1.0)
)
npi = Daedalus.DaedalusStructs.Npi([effect])

infection_re = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
infection_re.r0 = 5.0
@benchmark daedalus("Australia", infection_re, npi=npi, time_end=600.0)
```

```@example benchmarking_reactive_event
# shorter duration
infection_re2 = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
infection_re2.r0 = 5.0
@benchmark daedalus("Australia", infection_re2, npi=npi, time_end=100.0)
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

infection_te = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
infection_te.r0 = 3.0
@benchmark daedalus("Australia", infection_te, npi=timed_npi, time_end=600.0)
```

Timed events do not need $R_t$ logging and can be benchmarked without it.

```@example benchmarking_timed_event
infection_te2 = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
infection_te2.r0 = 3.0
@benchmark daedalus("Australia", infection_te2, npi=timed_npi, time_end=600.0, log_rt=false)
```
