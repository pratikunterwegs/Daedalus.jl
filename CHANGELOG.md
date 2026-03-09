# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
