# Plan: Data Infrastructure for Daedalus.jl

## Objective

Provide the same structured input data in Daedalus.jl as the R package `{daedalus.data}`, replacing
the current hardcoded Australia-specific values in `src/Data.jl` with a general, multi-country,
multi-pathogen data layer. This plan covers what datasets to include, where to store the raw files,
how to load and expose them, and in what order to implement them.

---

## Datasets in `{daedalus.data}`

| Dataset | R object | Contents | Dimensions |
|---------|----------|----------|------------|
| Country data | `country_data` | Demography, contact matrix, workers, GVA, hospital capacity | 197 countries |
| Economic contacts | `economic_contacts` | Within-sector contacts (`contacts_workplace`), cross-sector matrix (`contacts_between_sectors`) | 45-element vector, 45×45 matrix |
| Infection parameters | `infection_data` | Epi parameters for 7 pathogens | 7 pathogens × 12 parameters |
| Vaccination scenarios | `vaccination_scenario_data` | Start time, rate, uptake, efficacy, waning per scenario | 4 scenarios × 5 parameters |
| Closure / response strategies | `closure_data` | Per-sector openness coefficients per strategy | 4 strategies × 45 sectors |
| Sector names | `econ_sector_names` | Human-readable sector labels | 45 strings |
| Life value | `life_value` | Age-specific value of statistical life | 67 countries × 4 age groups |
| Life expectancy | `life_expectancy` | Age-specific life expectancy (years) | 67 countries × 4 age groups |
| Country GNI | `country_gni` | GNI per capita (international dollars, PPP) | 67 countries |
| Country names / codes | `country_names`, `country_codes_*` | ISO 3166 codes | 197 countries |

---

## Current State of `src/Data.jl`

- Hardcoded Australia demography, 4×4 contact matrix, worker counts.
- `worker_contacts()` uses synthetic data (explicitly flagged in its docstring).
- `consumer_worker_contacts()` returns a synthetic all-ones 45×4 matrix.
- `prepare_contacts()` bakes workplace and consumer-worker contacts into a single 49×49 matrix (to be separated per `plan_ode.md`).
- No multi-country or multi-pathogen support.

---

## Raw Data Files

The R package `{daedalus.data}` processes its data from CSV and Excel source files in
`inst/extdata/`. These are the canonical source files and should be **copied directly** into
this Julia package rather than regenerated.

### Files to bundle

| Source file | Format | Used for |
|-------------|--------|----------|
| `country_data.csv` | CSV | Country demography, contact matrices (16×16 age bins), worker counts |
| `sectorcontacts.csv` | CSV | Within-sector contact rates (`n_cnt` column per sector) |
| `sevenpathogens.csv` | CSV | All infection parameters for all 7 pathogens |
| `hospital_capacity.csv` | CSV | Spare hospital capacity by country (pre-processed) |
| `economic_closures.xlsx` → CSV | CSV (pre-export) | Sector openness coefficients by response strategy |
| `sector_gva_data.xlsx` → CSV | CSV (pre-export) | Gross value added per sector per country |

The two Excel files (`economic_closures.xlsx`, `sector_gva_data.xlsx`) should be exported to
CSV once using R or Python, then stored as CSV in the Julia package. This eliminates the xlsx
dependency from Julia.

### Destination

Create a new directory `src/data/` to hold these CSV files. This keeps data co-located with the
source module that loads it and avoids a separate `artifacts` infrastructure for data of this
scale.

```
src/
  data/
    country_data.csv        (from daedalus.data/inst/extdata/)
    sectorcontacts.csv      (from daedalus.data/inst/extdata/)
    sevenpathogens.csv      (from daedalus.data/inst/extdata/)
    hospital_capacity.csv   (from daedalus.data/inst/extdata/)
    economic_closures.csv   (pre-exported from economic_closures.xlsx)
    sector_gva_data.csv     (pre-exported from sector_gva_data.xlsx)
```

---

## Data Processing Logic

The R `data-raw/` scripts document how each dataset is derived. The Julia implementation
must replicate this logic at load time (or as a one-off preprocessing step). Key
transformations:

### Country data

**Demography** — `country_data.csv` has 21 five-year age bins (`Npop1`–`Npop21`).
Aggregate to 4 DAEDALUS groups:

| DAEDALUS group | Bins |
|----------------|------|
| `0-4` | 0–4 (bin 1) |
| `5-19` | 5–9, 10–14, 15–19 (bins 2–4) |
| `20-64` | 20–24, …, 60–64 (bins 5–13) |
| `65+` | 65–69, …, 100+ (bins 14–21) |

**Contact matrix** — `country_data.csv` has a 16×16 contact matrix (`CMaa`–`CMpp`).
Aggregate to 4×4 using population-weighted sums matching the demography aggregation.

**Workers** — columns `NNs1`–`NNs45` give workforce counts per sector. Fill missing values
with zero.

**GVA** — from `sector_gva_data.csv` (Daily sheet). Values are in million dollars per day.
Take absolute values (one country has a small negative typo).

**Hospital capacity** — from `hospital_capacity.csv`. Spare capacity =
`capacity_per_1000 × total_population × (1 − BOR/100)`.

### Economic contacts

`sectorcontacts.csv` column `n_cnt` gives 45 raw within-sector contact rates.
This becomes `contacts_workplace`. No normalisation by workforce count at this stage —
that division happens in the model layer (currently `worker_contacts()` in `Data.jl`).

`contacts_between_sectors` is a 45×45 zero matrix (diagonal zero, off-diagonal zero;
effectively no cross-sector contacts in the current model).

### Infection data

`sevenpathogens.csv` has rows for each parameter code and columns for each pathogen.
Key transformations from the R script:

- `sigma = 1 / Tlat`
- `gamma_Is = 1 / Tsr`
- `gamma_Ia = 1 / Tay`
- `gamma_H_recovery = 1 / Threc`
- `gamma_H_death = 1 / Thd`
- `rho = 1 / Ti`
- `eta = (IHR / p_sigma) / Tsh` (age-specific, aggregated to 4 groups by mean)
- `hfr = IFR / IHR` (age-specific)

IHR and IFR are given in 5-year age bins (`ihr1`–`ihr21`, `ifr1`–`ifr21`) and must be
averaged to the same 4-group structure as demography.

### Vaccination scenarios

Small enough to be hardcoded directly in Julia (the R script itself just defines
named vectors). No CSV required.

| Scenario | `start_time` (days) | `rate` (%pop/day) | `uptake_limit` (%) | `efficacy` (%) | `waning_period` (days) |
|----------|---------------------|-------------------|--------------------|----------------|------------------------|
| `none` | 365 | 1/7 | 40 | 50 | 270 |
| `low` | 300 | 2/7 | 50 | 50 | 270 |
| `medium` | 200 | 3/7 | 60 | 50 | 270 |
| `high` | 100 | 3.5/7 | 80 | 50 | 270 |

### Closure / response strategies

`economic_closures.csv` (exported from the `configurations` sheet of
`economic_closures.xlsx`) contains 45 rows (one per sector) and columns for each
configuration. The R script selects:

| Strategy | Implementation level used | Column in xlsx |
|----------|--------------------------|----------------|
| `elimination` | heavy | `lockdown` (strict version) |
| `economic_closures` | light | `economic_closures_(light)` |
| `school_closures` | light | `lockdown_school_closures_(light)` |
| `none` | — | All 1.0 (hardcoded) |

---

## Julia Module Design

### New module: `src/DataLoader.jl`

A new module `DataLoader` handles all dataset loading and caching. The existing `Data.jl`
retains Australia-specific convenience functions (backward compatibility) but delegates to
`DataLoader` where possible.

#### Struct definitions

```julia
module DataLoader

struct CountryData
    demography::Vector{Float64}        # 4 age groups
    contact_matrix::Matrix{Float64}    # 4×4
    workers::Vector{Int}               # 45 sectors
    gva::Vector{Float64}               # 45 sectors, million USD/day
    hospital_capacity::Float64         # spare beds
end

struct InfectionData
    r0::Float64
    sigma::Float64
    p_sigma::Float64
    epsilon::Float64
    gamma_Is::Float64
    gamma_Ia::Float64
    gamma_H_recovery::Float64
    gamma_H_death::Float64
    rho::Float64
    eta::Vector{Float64}               # 4 age groups
    hfr::Vector{Float64}              # 4 age groups
    ifr::Vector{Float64}              # 4 age groups
end

struct VaccinationScenario
    start_time::Float64
    rate::Float64
    uptake_limit::Float64
    efficacy::Float64
    waning_period::Float64
end

struct EconomicContacts
    contacts_workplace::Vector{Float64}             # 45 sectors
    contacts_between_sectors::Matrix{Float64}       # 45×45
end
```

#### Lazy-loaded module-level caches

Use `Ref{Union{Nothing, T}}` for lazy loading — parse CSV files on first access, cache
thereafter. This avoids loading all data at `using Daedalus` if only a subset is needed.

```julia
const _country_data_cache = Ref{Union{Nothing, Dict{String, CountryData}}}(nothing)
const _infection_data_cache = Ref{Union{Nothing, Dict{String, InfectionData}}}(nothing)
# etc.
```

#### Public API

```julia
# Country lookups
get_country(name::String)::CountryData
list_countries()::Vector{String}

# Pathogen lookups
get_pathogen(name::String)::InfectionData
list_pathogens()::Vector{String}

# Economic contacts (single global dataset)
get_economic_contacts()::EconomicContacts
get_sector_names()::Vector{String}

# Vaccination scenarios
get_vaccination_scenario(name::String)::VaccinationScenario
list_vaccination_scenarios()::Vector{String}

# Closure / response strategies
get_closure_strategy(name::String)::Vector{Float64}   # 45-element openness vector
list_closure_strategies()::Vector{String}
```

### Updates to `src/Data.jl`

Keep the existing Australia-specific functions (`australia_demography`,
`australia_contacts`, `australia_initial_state`, etc.) for backward compatibility, but:

1. Replace the synthetic `worker_contacts()` data with values from `sectorcontacts.csv`
   via `DataLoader.get_economic_contacts().contacts_workplace`.
2. Replace the all-ones `consumer_worker_contacts()` with appropriate data once
   real consumer-worker contact data is available (see Open Questions below).
3. Add a `prepare_community_contacts()` function as per `plan_ode.md`.

### Updates to `src/Daedalus.jl`

Include the new module:
```julia
include("DataLoader.jl")
```

---

## Implementation Priority

### Phase 1 — Unblock `plan_ode.md` (workplace contacts)

1. Copy `sectorcontacts.csv` to `src/data/`.
2. In `DataLoader.jl`, implement `get_economic_contacts()` to read
   `contacts_workplace` from the CSV (no struct needed yet; return a plain vector).
3. Update `Data.jl:worker_contacts()` to use real data from `sectorcontacts.csv`
   instead of the synthetic values.

### Phase 2 — Multi-pathogen support

4. Copy `sevenpathogens.csv` to `src/data/`.
5. Implement `InfectionData` struct and `get_pathogen()` / `list_pathogens()`.
6. Update `daedalus()` in `Model.jl` to accept an `InfectionData` argument
   (or keyword arguments matching its fields) as an alternative to manually
   specifying `r0`, `sigma`, etc.

### Phase 3 — Multi-country support

7. Copy `country_data.csv`, `hospital_capacity.csv`, `sector_gva_data.csv` to `src/data/`.
8. Implement demography and contact matrix aggregation from 21 bins → 4 groups and
   16×16 → 4×4.
9. Implement `CountryData` struct and `get_country()` / `list_countries()`.
10. Update `daedalus()` to accept a `CountryData` argument.
11. Refactor `australia_initial_state()` to call `get_country("Australia")` internally.

### Phase 4 — Response and vaccination data

12. Export `economic_closures.xlsx` → `economic_closures.csv` (one-time preprocessing).
13. Implement `get_closure_strategy()` for the 4 named response strategies.
14. Implement `VaccinationScenario` struct (hardcoded values; no CSV needed).
15. Update NPI handling in `Model.jl` / `Events.jl` to apply sector-specific openness
    vectors (links to `plan_ode.md` openness parameter).

### Phase 5 — Economic burden data (low priority)

16. Copy `life_expectancy.csv` and `value_life_year.csv`.
17. Implement life expectancy and GNI aggregation matching the R script logic.
18. Expose via `get_life_value()` and `get_life_expectancy()`.

---

## Dependencies to Add

| Package | Purpose |
|---------|---------|
| `CSV.jl` | Reading CSV files |
| `DataFrames.jl` | Tabular manipulation during data loading |

Both are standard Julia packages. Add to `Project.toml` under `[deps]`. `DataFrames.jl`
is only needed during loading/parsing; consider using plain matrix operations instead if
the dependency is undesirable.

---

## File Layout After Implementation

```
src/
  Daedalus.jl          (add DataLoader include)
  Constants.jl
  Data.jl              (keep; update worker_contacts + add prepare_community_contacts)
  DataLoader.jl        (new: structs + loading functions)
  DaedalusStructs.jl
  Helpers.jl
  Ode.jl
  Events.jl
  Model.jl
  Outputs.jl
  data/
    sectorcontacts.csv
    sevenpathogens.csv
    country_data.csv
    hospital_capacity.csv
    economic_closures.csv  (pre-exported from xlsx)
    sector_gva_data.csv    (pre-exported from xlsx)
```

---

## Open Questions

1. **Consumer-worker contacts**: The R package uses `contacts_consumer_worker`, which is
   in the `sec_contact_dist_UK.csv` file (age distribution of contacts across sectors for
   the UK). This file appears to be UK-specific and is not part of the `economic_contacts`
   dataset. The Julia package currently uses an all-ones placeholder. Should this remain
   synthetic, be loaded from `sec_contact_dist_UK.csv` as a fallback, or be derived from
   country demography data when available?

2. **Contact matrix aggregation**: The R aggregation of the 16×16 contact matrix to 4×4
   uses a population-weighted sum. The Julia implementation needs to reproduce this exactly
   to match R outputs. Should this be validated against known Australia values as a test?

3. **GVA data**: The `sector_gva_data.xlsx` file is not currently used in the Julia ODE model
   (only epi dynamics, not economic outcomes, are modelled). Include it in Phase 3 only as
   stored data, or defer entirely until economic analysis is needed?

4. **Data bundling strategy**: For a small-to-medium package, bundling CSV files directly
   in `src/data/` is simplest. If the package is ever published and the data files grow
   large (e.g. all 197 countries' contact matrices), consider migrating to Julia's
   `Artifacts.jl` system for versioned, lazily-downloaded data.

5. **Test coverage**: Should a data-validation test suite mirror the `checkmate` assertions
   in the R data-raw scripts (e.g. verify that each country's contact matrix is 4×4,
   non-negative, finite)? This would guard against CSV corruption.
