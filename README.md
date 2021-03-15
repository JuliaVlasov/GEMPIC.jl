# GEMPIC.jl

Geometric ElectroMagnetic Particle-In-Cell Methods

**NOTE: This package is still very much under development.**


<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) 
-->
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
[![Build Status](https://travis-ci.org/JuliaVlasov/GEMPIC.jl.svg?branch=master)](https://travis-ci.org/JuliaVlasov/GEMPIC.jl)
[![codecov](https://codecov.io/gh/JuliaVlasov/GEMPIC.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaVlasov/GEMPIC.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliavlasov.github.io/GEMPIC.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliavlasov.github.io/GEMPIC.jl/dev)

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

This is a translation from a Fortran code written by :

- Yaman Güçlü
- Katharina Kormann  
- Benedikt Perse
- Eric Sonnendrücker
- Edouardo Zoni

from [Max-Planck-Institut fur Plasmaphysik - Garching (Germany)](https://www.ipp.mpg.de/4098496/kgkm)
