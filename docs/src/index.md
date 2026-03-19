```@meta
CurrentModule = Daedalus
```

# Daedalus.jl

This is some minimal documentation for [Daedalus.jl](https://github.com/jameel-institute/Daedalus.jl).

[![License:MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jameel-institute.github.io/Daedalus.jl/dev/)
[![Build Status](https://github.com/jameel-institute/Daedalus.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jameel-institute/Daedalus.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jameel-institute/Daedalus.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jameel-institute/Daedalus.jl)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

_Daedalus.jl_ is a Julia package that aims to mirror the [R package {daedalus}](https://github.com/jameel-institute/daedalus).

This section shows how to run the Daedalus Julia model.
Functionality compared to the R package is limited.

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

## Running Daedalus for a specific country

Pass a country name string directly to `daedalus` to use country-specific demography, contact patterns, and workforce data.

```@example uk_daedalus
using Daedalus
using Plots

# Run the model using UK demography and contact patterns
infection_uk = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
infection_uk.r0 = 2.5
data_uk = daedalus("United Kingdom", infection_uk, time_end=600.0)

times_uk = Daedalus.Outputs.get_times(data_uk)
exposed_uk = Daedalus.Outputs.get_values(data_uk, "E", 1)
hosp_uk = Daedalus.Outputs.get_values(data_uk, "H", 1)

plot(times_uk, exposed_uk, label = "exposed")
plot!(times_uk, hosp_uk, label = "hospitalised")
xlabel!("Time (days)")
ylabel!("# individuals")
title!("United Kingdom — SEIR dynamics")
```

```@example uk_daedalus
# Effective reproduction number over time
rt_uk = Daedalus.Outputs.get_values(data_uk, "Rt", 1)
plot(times_uk, rt_uk, label = "Rt", color = :red)
hline!([1.0], linestyle = :dash, color = :black, label = "Rt = 1")
xlabel!("Time (days)")
ylabel!("Rt")
title!("United Kingdom — effective reproduction number")
```

See [Country and pathogen data](@ref) for a full overview of the bundled
data and how to use pathogen-specific parameters.
