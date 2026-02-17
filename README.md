# Daedalus

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pratikunterwegs.github.io/Daedalus.jl/dev/)
[![Build Status](https://github.com/pratikunterwegs/Daedalus.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pratikunterwegs/Daedalus.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/pratikunterwegs/Daedalus.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/pratikunterwegs/Daedalus.jl)

_Daedalus.jl_ is a Julia package that aims to mirror the [R package {daedalus}](https://github.com/jameel-institute/daedalus).

**Note** that this is a personal project, and comes with no current or future support.
This documentation section is intended as a learning experience (for me) in writing Julia package documentation.

## Installation

_Daedalus.jl_ can be installed from GitHub using the Julia package manager _Pkg.jl_.

```julia
using Pkg
Pkg.add(url="git@github.com:pratikunterwegs/Daedalus.jl.git")
```

## Quick start

```julia
using Daedalus

daedalus(r0=5.0, time_end=600.0)
```

## Related projects

- The Daedalus integrated epi-econ model provided by the R package [_daedalus_](https://jameel-institute.github.io/daedalus/).

## Help

To report a bug, request a feature, or just start a discussion, [please open an issue](https://github.com/pratikunterwegs/Daedalus.jl/issues/new).
