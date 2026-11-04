# Strong Landau damping in 2d2v

The 2d2v example evolves two spatial coordinates and two velocity components
with periodic boundaries. It uses the existing 2D Maxwell solver and the
Hamiltonian Strang splitting. This is the closed invariant subsystem with
electric components ``E_x,E_y`` and magnetic component ``B_z``. The other
Maxwell components ``E_z,B_x,B_y`` are initialized to zero and remain zero.

The initial spatial density uses the same amplitude and wave number as the
repository's strong 1d2v reference case:

```math
n_0(x,y)=1+0.5\cos(x/2)+0.5\cos(y/2),
```

on ``[0,4\pi]^2``, with a unit-temperature centered Maxwellian in both velocity components. The
anisotropic spline spaces match the 2D Maxwell solver: ``E_x`` uses degrees
``(p-1,p)``, ``E_y`` uses ``(p,p-1)``, and ``B_z`` uses ``(p-1,p-1)``. Charge
density uses degree ``(p,p)`` with Galerkin scaling, as required for the weak
charge right-hand side consumed by the Poisson solver. The example plots the
logarithm of the total electric energy ``E_x+E_y``; the strong nonlinear
evolution need not be a straight exponential or decrease monotonically at
every time.

Run the example from the repository root:

```sh
julia --project examples/strong_landau_damping_2d2v.jl
```

The `run_strong_landau_damping_2d2v` function accepts `steps`, `n_particles`,
`nx`, `ny`, and `dt` keyword arguments and returns a DataFrame of particle and
field diagnostics.
