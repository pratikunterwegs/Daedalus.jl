# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.0.10] - 20266-04-09

### Changed (Breaking)
- Simplified public API by hiding internal trigger state representation:
  - Removed from exports: `ReactiveEffect`, `TimedEffect`, `StateData`, `CtStateData`, `DsStateData`, `get_indices`
  - Added to exports: `Effect` (abstract type), `ParamEffect` (concrete effect type), `Trigger` (abstract type), `ReactiveTrigger`, `TimeTrigger` (concrete trigger types)
  - Internal state tracking now uses abstract `Trigger` type instead of concrete `StateData` subtypes
- Trigger constructor parameter order changed (breaking change for direct instantiation):
  - `ReactiveTrigger(name, value)` → `ReactiveTrigger(value, name)` (value now first positional argument)
  - `TimeTrigger(name="time", value)` → `TimeTrigger(value, name="time")` (value now first positional argument)
- Removed convenience constructor: `ReactiveEffect(target, func, reset_func; on=(...), off=(...))` — use full `ParamEffect` constructor with explicit `ReactiveTrigger`/`TimeTrigger` objects instead
- Removed public function: `get_indices(x::StateData)` — internal helper no longer exported
- Refactored `make_param_changer(eff::Effect)` and `make_param_reset(eff::Effect)` to dispatch on `Trigger` type using `isa(eff.trigger_on, TimeTrigger)` instead of string comparison on trigger names

## [0.0.9] - 2026-04-01

### Changed (Breaking)
- Unified NPI types: `Npi` is now the single container for all parameter effects (both reactive and timed)
  - `abstract type Effect` renamed to `abstract type ParamEffect`
  - `ParamEffect` struct renamed to `ReactiveEffect` — state-dependent interventions with independent on/off triggers
  - New `TimedEffect <: ParamEffect` struct for time-limited interventions (replaces `TimedNpi`)
    - Constructor: `TimedEffect(target, func, reset_func, start_time, end_time)`
    - Example: `TimedEffect(:beta, x -> x .* 0.7, x -> x ./ 0.7, 10.0, 40.0)`
  - `Npi` now accepts both `ReactiveEffect` and `TimedEffect` objects via `Npi(effects::AbstractVector{<:ParamEffect})`
    - Unified callback handling: `make_events(npi, savepoints)` dispatches per effect type
  - `TimedNpi` struct removed; use `Npi([TimedEffect(...)])` instead
  - Removed public functions: `n_phases`, `total_duration`, `make_timed_npi_callbacks`
  - Removed public type: `TimedNpi`
- Behavior changes:
  - `result.saves` now only contains `SavedValues` for `ReactiveEffect` entries (skips `TimedEffect`)
  - `Union{Npi, TimedNpi, Nothing}` type annotations throughout codebase replaced with `Union{Npi, Nothing}`
  - `make_param_changer(eff::ParamEffect)` and `make_param_reset(eff::ParamEffect)` now dispatch on specific effect types
  - `make_save_events(npi::Npi)` skips `TimedEffect` entries (no state-based saving needed for timed interventions)

## [0.0.8] - 2026-03-31

### Changed (Breaking)
- Added struct `ParamEffect` to define NPI target parameter and launch and end triggers
  - New constructor: `ParamEffect(target, func, reset_func; on=(...), off=(...))`
  - `func`: Transformation function applied when effect is active (e.g., `x -> x .* 0.4` for 60% reduction)
  - `reset_func`: Inverse function to restore parameter when effect deactivates (e.g., `x -> x ./ 0.4`)
  - Example: `ParamEffect(:beta, x -> x .* 0.4, x -> x ./ 0.4; on=("H", 5000.0), off=("Rt", 1.0))`
- Struct `Npi` simplified to a container of `Vector{ParamEffect}`
  - Removed fields: `params` (NpiData), `saved_values`, `ison` — all now per-effect
  - New primary constructor: `Npi(effects::Vector{ParamEffect})` for custom per-effect triggers
- `make_events(npi::Npi, savepoints)`: old signature `make_events(npi, on_func, off_func, savepoints)` removed
- `make_save_events(npi::Npi, savepoints)`" now returns `Vector{SavingCallback}` (one per effect) instead of single callback
- Internal logic for parameter modification moved from `make_param_changer(npi)` into `make_events`
- `result.saves` is now a `Vector{SavedValues}` (one per effect) instead of a single SavedValues
- Function `get_coef`: Deprecated and now throws an error; marked for removal
- Flexible `Npi` parameter specifications: `Npi` now accepts flexible parameter names and transform functions via new constructors:
  - Per-parameter functions (Vector of Pairs): `Npi(threshold::Float64, [:param1 => func1, :param2 => func2])` — different transformation function for each parameter
  - Per-parameter functions (Dict): `Npi(threshold::Float64, Dict(:param1 => func1, :param2 => func2))` — Dict-based variant of per-parameter functions
- New `make_param_changer(npi::Npi)` dispatch: Replaces hardcoded parameter modification logic. Reads parameter specifications from the `Npi` struct and applies all effects in a single callback.
- New `make_param_reset(npi::Npi)` dispatch: Companion to `make_param_changer` that resets all modified parameters to their original values.

### Changed
- `Npi` struct internals: Replaced `coefs::NamedTuple` field with `effects::Vector{ParamEffect}` to support flexible parameter specifications.
- Removed `get_coef` from public API: No longer exported; `make_param_changer(npi::Npi)` and `make_param_reset(npi::Npi)` handle coefficient extraction internally.


## [0.0.7] - 2026-03-19

- Transfer repository to Jameel Institute @jameel-institute oragnisation
- Added small script `docs/sync_readme.jl` to update package version in website index and Readme.md

## [0.0.6] - 2026-03-16

### Changed (Breaking)
- **`daedalus()` function signature refactored**: second positional argument now takes infection parameters (`infection`) instead of scalar/vector `r0`. Users must pass parameters via `InfectionData` object rather than individual keyword arguments
  - **Old interface**: `daedalus(country, r0::Float64; sigma=..., epsilon=..., eta=..., ...)`
  - **New interface**: `daedalus(country, infection; npi=..., log_rt=..., time_end=..., increment=..., n_threads=...)`
- Removed all infection-parameter keyword arguments: `sigma`, `p_sigma`, `epsilon`, `rho`, `eta`, `omega`, `gamma_Ia`, `gamma_Is`, `gamma_H`, `nu`
  - All epidemiological parameters now encapsulated in `InfectionData` object
  - Users customize parameters by fetching `InfectionData` and modifying fields before calling `daedalus()`

### Added
- New dispatch methods for `daedalus()`:
  1. String pathogen name: `daedalus(country, "sars-cov-2 delta"; ...)`
  2. Single `InfectionData`: `daedalus(country, infection_obj; ...)`
  3. Vector `InfectionData`: `daedalus(country, [inf1, inf2, ...]; ...)`
- `extract_infection_params()` helper function to extract and expand epidemiological parameters from `InfectionData`
- `InfectionData` is mutable, allowing users to customize parameters post-fetch: `inf = get_pathogen("sars-cov-2 delta"); inf.r0 = 2.5`
- **Infection names are normalized to lowercase**: all pathogen names are stored and looked up as lowercase strings (e.g., `"sars-cov-2 delta"`, `"influenza 2009"`)

### Migration Guide
```julia
# Old (no longer works):
result = daedalus("Australia", 2.5, sigma=0.217, epsilon=0.58, time_end=200.0)

# New (string pathogen, lowercase names):
result = daedalus("Australia", "sars-cov-2 delta", time_end=200.0)

# New (custom infection):
infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
infection.r0 = 2.5
result = daedalus("Australia", infection, time_end=200.0)

# New (vector of infections):
infections = [
    Daedalus.DataLoader.get_pathogen("sars-cov-2 delta"),
    Daedalus.DataLoader.get_pathogen("influenza 2009")
]
results = daedalus("Australia", infections, time_end=200.0)
```

## [0.0.5] - 2026-03-09

### Changed
- `daedalus()` function signature: `country` is now the first positional argument and `r0` is the second positional argument, enabling method dispatch on `r0` type (scalar `Float64` vs vector `Vector{Float64}`) across two separate implementations in `src/Model.jl` and `src/Ensemble.jl`
- All documentation examples and benchmarks updated to reflect the new calling convention: `daedalus(country, r0, ...)` instead of previous keywords-first approach
- All test files updated to use new positional argument signature for daedalus calls
- Function calls now use positional arguments: `daedalus("Australia", 2.5, time_end=200.0)` instead of previous keyword-based calling conventions
- Function `get_ngm()` requires transmission rate beta and not the $R_0$; function `get_beta()` is now vectorised over multiple values of $R_0$. The use case is generating multiple NGMs for ensemble runs without running `get_beta()` an equal number of times. `get_ngm()` has a method for a vector of `beta()`.

### Added
- Multiple dispatch implementation for `daedalus()`: scalar and vector R0 inputs are now handled via distinct function methods
- Vector R0 dispatch in file `src/Ensemble.jl`: `daedalus(country, r0::Vector{Float64}; ...)` runs multiple R0 values in a single call
- Implementation of `SciMLBase.EnsembleProblem` in `src/Model.jl` (`daedalus_internal` function): uses `EnsembleThreads()` solver with `prob_func` wrapper to efficiently orchestrate multi-run ODE solving. The ensemble approach reuses a base ODE problem and remakes it for each trajectory with its corresponding parameters, enabling automatic thread-safe parallel execution across multiple r0 values without explicit locking
- Helper functions `prepare_shared_data()` and `daedalus_internal()` exported from Model.jl for use by ensemble dispatch
- Comprehensive function documentation improvements: added or enhanced docstrings for all exported functions in `Helpers.jl`, `Data.jl`, `Events.jl`, `Model.jl`, and `Ode.jl` with argument lists and return type annotations

## [0.0.4] - 2026-03-02

### Added
- Documentation page `docs/src/settings.md` explaining the multiple-contact-settings feature: how to assign a `Vector{Matrix{Float64}}` to `CountryData.contact_matrix`, how `contacts3d` stacks them into a 3D array, and how beta calibration uses `total_contacts` (sum of all matrices)
- Tests for multiple contact settings: `get_settings` count, `contacts3d` shape, `total_contacts` element-wise sum, model execution with two settings, and calibration-equivalence check (two equal settings produces the same epidemic as one setting for the same R0)

### Fixed
- Lowered `Statistics` and `LinearAlgebra` compat bounds to `1.10.0` (matching Julia 1.10 LTS stdlib versions) so the package resolves correctly on Julia 1.10 LTS

### Changed

- `CountryData` struct now accepts a `Vector` of contact matrices as `Matrix{Float64}` for multiple contact settings. Helper functions process this list, or a single `Matrix`, to give total contacts where needed including `Helpers.get_beta` and `Helpers.get_ngm`.
- Moved away from using `StaticArrays` for contact matrices as operating on them was slower than using regular arrays.

### Added
- Docstring for `Helpers.weighted_slice_sum!` explaining the tensor contraction algorithm, arguments, and performance notes
- Docstring for `Data.total_contacts` explaining the dispatch on single vs. vector-of-matrices input
- Expanded docstring for `Data.contacts3d` explaining the 3D stacking, the `K=1` reshape fallback, and the role of the third dimension in the ODE
- Tests for `Helpers.weighted_slice_sum!` covering unit weights, slice selection, scalar scaling, zero weights, in-place overwrite, and agreement with a reference loop

## [0.0.3] - 2026-02-27

### Changed
- `daedalus` now accepts `country` as either a `String` or a `DataLoader.CountryData` struct; a `String` is resolved to `CountryData` via `DataLoader.get_country` at the start of the function, making both call styles equivalent
- Updated `test/test_basic.jl` and `test/test_eigenvalue.jl` to replace removed `Data.australia_contacts()` calls with `DataLoader.get_country("Australia").contact_matrix`; replaced zero-arg `Data.prepare_contacts()` call with `Data.prepare_contacts("Australia")`

## [0.0.2] - 2026-02-25

### Fixed
- `prepare_demog(cd::CountryData)` now clamps worker counts to a minimum of 1,
  preventing division-by-zero when the result is used as a scaling denominator.
  30 countries (Australia, Belgium, Brunei, Cambodia, Chile, China, Costa Rica,
  Cyprus, Estonia, Finland, Hong Kong, Iceland, Japan, Kazakhstan, Laos, Latvia,
  Luxembourg, Malaysia, Malta, Mexico, Morocco, Myanmar, New Zealand, Portugal,
  Romania, Rwanda, Singapore, Slovenia, Switzerland, Tunisia) had at least one
  sector with zero workers in the data, causing `Inf` in the scaled contact
  matrix, which propagated to `NaN` in the ODE step-size calculation and
  immediate solver exit with `dt_NaN` warnings.

### Changed
- `daedalus` now requires a `country` string directly instead of separate `initial_state`, `contacts`, and `cw` arguments; all tests, examples, and documentation updated accordingly
- Simplified UK example in `docs/src/index.md` to use `daedalus(country="United Kingdom", ...)`
- Updated `docs/src/country_data.md` to reflect that `daedalus` accepts `country` directly

### Added
- Documentation page for country and pathogen data (`docs/src/country_data.md`)
- `DataLoader` module added to function reference autodocs
- Basic model and helper tests extracted into `test/test_basic.jl`

### Fixed
- Added explicit `du[end] = 0.0` in `daedalus_ode!` to prevent undefined Rt derivative between callback updates
- Changed `Ia * p.epsilon` to `Ia .* p.epsilon` in ODE for consistency with broadcasting conventions
- Removed `cm_scaling .* p.contacts` StaticArrays broadcasting failure; replaced with `sum(p.contacts, dims=3)[:,:,1]`

## [0.0.1] - 2026-02-23

### Added
- Core epidemiological model mirroring the R package {daedalus}
- ODE-based compartmental disease transmission model with demographic groups
- Effective reproduction number (Rt) calculation and logging at specified timesteps
- Next Generation Matrix (NGM) method for Rt calculation
- Power iteration for faster Rt calculation
- Flexible event handling system with both timed and reactive (state-dependent) events
- Non-pharmaceutical intervention (NPI) modeling with timed and reactive triggers
- Time dependent NPIs as a struct `TimedNpi`
- Contact matrix support for modeling population mixing patterns
- Option to toggle contact matrix scaling
- Functions for dynamic parameter modification and reset during simulations
- Support for both increasing and decreasing threshold detection in event callbacks
- Worker contact modeling with static vector optimization
- Documentation with basic usage examples and plots
- Documentation workflow via GitHub Actions
- Data layer (`DataLoader`) with lazy-loaded country, pathogen, economic contacts, closure strategy, and vaccination scenario data mirroring the R {daedalus} data package
- Bundled CSV data files: country demographics, hospital capacity, sector GVA, sector contacts, economic closure strategies, and seven pathogen parameter sets
- Output timeseries access function for structured post-simulation analysis

### Changed
- Simplified NPI handling and data structures
- Improved event logic for better state-dependent triggering
- Unified contact calculation using single contact matrix approach
- Optimized dependency management for lighter package footprint
- Updated ODE system to account for worker-specific transmission dynamics
- Refined beta (transmission rate) calculation for arbitrary contact matrix sizes
- Enhanced model interface to return NPIs and remove unused inputs
- Applied JuliaFormatter across codebase for consistent style

### Fixed
- Corrected Rt calculation method
- Fixed ODE formulation for accurate disease dynamics
- Resolved issues with `SavedValues` type handling in callbacks
- Fixed error when output timebin is an exact factor of tmax

### Technical Details
- Package version: 0.0.1
- Julia compatibility: 1.6.7+
- Key dependencies: OrdinaryDiffEq.jl, DiffEqCallbacks.jl, StaticArrays.jl, CSV.jl, DataFrames.jl
