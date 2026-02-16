# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Core epidemiological model mirroring the R package {daedalus}
- ODE-based compartmental disease transmission model with demographic groups
- Effective reproduction number (Rt) calculation and logging at specified timesteps
- Next Generation Matrix (NGM) method for Rt calculation
- Flexible event handling system with both timed and reactive (state-dependent) events
- Non-pharmaceutical intervention (NPI) modeling with timed and reactive triggers
- Contact matrix support for modeling population mixing patterns
- Option to toggle contact matrix scaling
- Functions for dynamic parameter modification and reset during simulations
- Support for both increasing and decreasing threshold detection in event callbacks
- Worker contact modeling with static vector optimization
- Documentation with basic usage examples and plots
- Time dependent NPIs as a struct `TimedNpi`

### Changed
- Simplified NPI handling and data structures
- Improved event logic for better state-dependent triggering
- Unified contact calculation using single contact matrix approach
- Optimized dependency management for lighter package footprint
- Updated ODE system to account for worker-specific transmission dynamics
- Refined beta (transmission rate) calculation for arbitrary contact matrix sizes
- Enhanced model interface to return NPIs and remove unused inputs

### Fixed
- Corrected Rt calculation method
- Fixed ODE formulation for accurate disease dynamics
- Resolved issues with `SavedValues` type handling in callbacks

### Technical Details
- Package version: 1.0.0-DEV
- Julia compatibility: 1.6.7+
- Key dependencies: OrdinaryDiffEq.jl, DiffEqCallbacks.jl, StaticArrays.jl
- Development status: Work in Progress (WIP)
