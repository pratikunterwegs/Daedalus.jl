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

delta = Daedalus.DataLoader.get_pathogen("SARS-CoV-2 delta")

delta
```

## Preparing country inputs for the model

Pass a country name string directly to [`daedalus`](@ref) to run the model with country-specific demography, contacts, and workforce data:

```@example country_inputs
using Daedalus

result = daedalus(country="United Kingdom", r0=2.5, time_end=600.0)
```

You can also pass a [`DataLoader.CountryData`](@ref) struct directly.
This is useful when you want to pre-fetch or modify country data before running the model:

```@example country_inputs
# Pre-fetch country data and pass the struct directly
uk = Daedalus.DataLoader.get_country("United Kingdom")
result2 = daedalus(country=uk, r0=2.5, time_end=600.0)
```

The `Data` sub-module also exposes lower-level functions for inspecting or manipulating the country arrays directly.
Each function accepts either a country name `String` or a `CountryData` struct:

| Function | Returns |
|---|---|
| `initial_state(country)` | Initial epidemic state array (N_GROUPS × N_COMPARTMENTS × N_VAX_STRATA) |
| `prepare_contacts(country)` | Scaled 49×49 contact matrix |
| `worker_contacts(country)` | Per-capita within-sector contact rates (length-45 vector) |
| `prepare_demog(country)` | 49-element population vector |
