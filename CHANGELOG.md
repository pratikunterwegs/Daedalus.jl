# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- `daedalus` now accepts `country` as either a `String` or a `DataLoader.CountryData` struct; a `String` is resolved to `CountryData` via `DataLoader.get_country` at the start of the function, making both call styles equivalent

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
