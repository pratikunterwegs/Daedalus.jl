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

## Benchmarking

```@example benchmarking
# benchmark for a typical daedalus run of 600 days
using Daedalus
using BenchmarkTools

@benchmark daedalus()
```
