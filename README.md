# GEMPIC.jl

Geometric ElectroMagnetic Particle-In-Cell Methods

[![CI](https://github.com/JuliaVlasov/GEMPIC.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaVlasov/GEMPIC.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaVlasov/GEMPIC.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaVlasov/GEMPIC.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliavlasov.github.io/GEMPIC.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliavlasov.github.io/GEMPIC.jl/dev)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

A Julia implementation of the [GEMPIC](https://arxiv.org/abs/1609.03053)

## Installation

In a Julia session switch to `pkg>` mode to add `GEMPIC`:

```julia
julia>] # switch to pkg> mode
pkg> add GEMPIC
```

When finished, make sure that you're back to the Julian prompt (`julia>`)
and bring `GEMPIC` into scope:

```julia
julia> using GEMPIC
```

## Credits

This is a translation from Fortran code [selalib](https://github.com/selalib/selalib) written by :

- Yaman Güçlü
- Katharina Kormann  
- Benedikt Perse
- Eric Sonnendrücker
- Edouardo Zoni

from [Max-Planck-Institut fur Plasmaphysik - Garching (Germany)](https://www.ipp.mpg.de/4098496/kgkm)
