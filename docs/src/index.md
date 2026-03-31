```@meta
CurrentModule = Daedalus
```

# Daedalus.jl

[![License:MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/version-0.0.8-aquamarine.svg)](https://jameel-institute.github.io/Daedalus.jl/dev/)
[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jameel-institute.github.io/Daedalus.jl/dev/)
[![Build Status](https://github.com/jameel-institute/Daedalus.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jameel-institute/Daedalus.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jameel-institute/Daedalus.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jameel-institute/Daedalus.jl)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

_Daedalus.jl_ is a Julia package that aims to mirror the [R package {daedalus}](https://github.com/jameel-institute/daedalus).

_Daedalus.jl_ is developed at the [Jameel Institute](https://www.imperial.ac.uk/jameel-institute/) at Imperial College London as part of the [Jameel Institute-Kenneth C. Griffin Initiative for the Economics of Pandemic Preparedness (EPPI)](https://new.express.adobe.com/webpage/TXLBkz1sN9FI5?), in collaboration with the [RESIDE research software engineering team](https://reside-ic.github.io/about/).

**Note that** functionality compared to the R package is still limited.

## Installation

_Daedalus.jl_ can be installed from GitHub using the Julia package manager _Pkg.jl_.

```julia
using Pkg
Pkg.add(url="git@github.com:jameel-institute/Daedalus.jl.git")
```

## Quick start

```@example basic_daedalus
using Daedalus
using Plots

# pump up r0 to get peak within 50 days
infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
infection.r0 = 5.0
data = daedalus("Canada", infection, time_end=600.0);

# plot exposed group
times = Daedalus.Outputs.get_times(data)
# functioning get_values() bins into 90 days by default, set to 1 for no binning
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

## Related projects

- The Daedalus integrated epi-econ model provided by the R package [_daedalus_](https://jameel-institute.github.io/daedalus/).

## Help

To report a bug, request a feature, or just start a discussion, [please open an issue](https://github.com/jameel-institute/Daedalus.jl/issues/new).
