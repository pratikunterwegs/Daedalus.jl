```@meta
CurrentModule = Daedalus
```

# Daedalus.jl

This is some minimal documentation for [Daedalus.jl](https://github.com/pratikunterwegs/Daedalus.jl).

**Note** that this is a personal project, and comes with no current or future support.
This documentation section is intended as a learning experience (for me) in writing Julia package documentation.

[![License:MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pratikunterwegs.github.io/Daedalus.jl/dev/)
[![Build Status](https://github.com/pratikunterwegs/Daedalus.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pratikunterwegs/Daedalus.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/pratikunterwegs/Daedalus.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/pratikunterwegs/Daedalus.jl)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

_Daedalus.jl_ is a Julia package that aims to mirror the [R package {daedalus}](https://github.com/jameel-institute/daedalus).

This section shows how to run the Daedalus Julia model.
Functionality compared to the R package is limited.

```@example basic_daedalus
using Daedalus
using Plots

# pump up r0 to get peak within 50 days
data = daedalus(r0=5.0, time_end=600.0);

# plot exposed group
times = Daedalus.Outputs.get_times(data)
# functiong get_values() bins into 90 days by default, set to 1 for no binning
exposed = Daedalus.Outputs.get_values(data, "E", 1)
hosp = Daedalus.Outputs.get_values(data, "H", 1)

# plot the output to see lag in hospitalisations
plot(times, exposed, label="exposed")
plot!(times, hosp, label="hosp")
xlabel!("Time (days)")
ylabel!("# individuals")
```

```@example basic_daedalus
# plot recorded Rt
vecRt = Daedalus.Outputs.get_values(data, "Rt", 1)
plot(times, vecRt, label="Rt")
xlabel!("Time (days)")
ylabel!("Rt")
```
