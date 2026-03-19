```@meta
CurrentModule = Daedalus
```

# Country and pathogen data

_Daedalus.jl_ bundles demographic, economic, and epidemiological data that mirrors the R package [{daedalus.data}](https://github.com/jameel-institute/daedalus.data).
All data access is provided through the `DataLoader` sub-module, which lazy-loads CSV files on first access and caches results for subsequent calls.

## Listing available data

```@example data_listing
using Daedalus

# All countries with bundled demographic + economic data
countries = Daedalus.DataLoader.list_countries()
println("$(length(countries)) countries available")
println("First five: ", countries[1:5])
```

```@example data_listing
# All seven reference pathogens
Daedalus.DataLoader.list_pathogens()
```

```@example data_listing
# Pandemic response (economic closure) strategies
Daedalus.DataLoader.list_closure_strategies()
```

```@example data_listing
# Pre-pandemic vaccine investment scenarios
Daedalus.DataLoader.list_vaccination_scenarios()
```

## Country data

[`DataLoader.get_country`](@ref) returns a [`DataLoader.CountryData`](@ref) struct containing demography, a 4×4 social contact matrix, sector workforce counts, gross value added (GVA) per sector, and hospital capacity.

```@example country_struct
using Daedalus

uk = Daedalus.DataLoader.get_country("United Kingdom")

println("Demography (4 age groups): ", uk.demography)
println("Hospital capacity (spare beds): ", uk.hospital_capacity)
println("Number of economic sectors: ", length(uk.workers))
```

```@example country_struct
# 4×4 social contact matrix
uk.contact_matrix
```

## Pathogen data

[`DataLoader.get_pathogen`](@ref) returns an [`DataLoader.InfectionData`](@ref) struct with all epidemiological parameters derived from the seven reference pathogens.

```@example pathogen_struct
using Daedalus

delta = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")

delta
```

## Preparing country inputs for the model

Pass a country name string and a pathogen name to [`daedalus`](@ref) to run the model with country-specific demography, contacts, and workforce data:

```@example country_inputs
using Daedalus

result = daedalus("United Kingdom", "sars-cov-2 delta", time_end=600.0)
```

You can also pass a [`DataLoader.CountryData`](@ref) struct directly and customize infection parameters.
This is useful when you want to pre-fetch or modify country or infection data before running the model:

```@example country_inputs
# Pre-fetch country data and pass the struct directly
uk = Daedalus.DataLoader.get_country("United Kingdom")
infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
infection.r0 = 2.5
result2 = daedalus(uk, infection, time_end=600.0)
```

The `Data` sub-module also exposes lower-level functions for inspecting or manipulating the country arrays directly.
Each function accepts either a country name `String` or a `CountryData` struct:

| Function | Returns |
|---|---|
| `initial_state(country)` | Initial epidemic state array (N_GROUPS × N_COMPARTMENTS × N_VAX_STRATA) |
| `prepare_contacts(country)` | Scaled 49×49 contact matrix |
| `worker_contacts(country)` | Per-capita within-sector contact rates (length-45 vector) |
| `prepare_demog(country)` | 49-element population vector |

## Running Daedalus with country and infection structs

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
