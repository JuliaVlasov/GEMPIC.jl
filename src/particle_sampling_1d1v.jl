using Sobol, Random, Distributions

export ParticleSampler

"""
    ParticleSampler{D,V}( sampling_type, symmetric, dims, n_particles)

Particle initializer class with various functions to initialize a particle.

- `sampling_type` : `:random` or `:sobol`
- `symmetric` : `true` or `false`
- `n_particles` : number of particles
"""
struct ParticleSampler{D,V}
    sampling_type::Symbol
    dims::Tuple{Int,Int}
    n_particles::Int
    symmetric::Bool
    seed::Int

    function ParticleSampler{D,V}(
        sampling_type::Symbol, symmetric::Bool, n_particles::Int, seed::Int=1234
    ) where {D,V}
        if !(sampling_type in [:random, :sobol])
            throw(ArgumentError("Sampling type $sampling_type 
                      not implemented"))
        end

        dims = (D, V)

        # Make sure that the particle number is conforming 
        # with symmetric sampling (if necessary)
        if symmetric
            ncopies = 2^(sum(dims))
            np = mod(n_particles, ncopies)
            if np != 0
                n_particles += np
            end
        else
            ncopies = 1
        end

        return new(sampling_type, dims, n_particles, symmetric, seed)
    end

    function ParticleSampler{D,V}(
        sampling_type::Symbol, n_particles::Int, seed::Int=1234
    ) where {D,V}
        if !(sampling_type in [:random, :sobol])
            throw(ArgumentError("Sampling type $sampling_type 
                      not implemented"))
        end

        dims = (D, V)
        symmetric = false

        return new(sampling_type, dims, n_particles, symmetric, seed)
    end
end

export sample!

"""
     newton(r, α, k)

Function to solve ``P(x) - r = 0`` where ``r \\in [0, 2π/k]``

where ``P`` is the cdf of ``f(x) = 1 + α \\cos(k x)``
     
"""
function newton(r, α, k)
    x0, x1 = 0.0, 1.0
    r *= 2π / k
    while (abs(x1 - x0) > 1e-12)
        p = x0 + α * sin(k * x0) / k
        f = 1 + α * cos(k * x0)
        x0, x1 = x1, x0 - (p - r) / f
    end
    x1
end

function sample!( pg::ParticleGroup{1,1}, ps::ParticleSampler, df::AbstractCosGaussian, mesh::OneDGrid)

    α = df.params.α[1] 
    k = df.params.k[1][1]
    σ = df.params.σ[1][1]
    sample!( pg, α, k, σ, mesh)

end

"""
    sample!( pg::ParticleGroup{1,1}, α, k, σ, mesh::OneDGrid)

Sampling from a probability distribution to initialize a Landau damping in
1D1V space.

```math
f_0(x,v,t) = \\frac{n_0}{\\sqrt{2π} v_{th}} ( 1 + \\alpha cos(k_x x)) exp( - \\frac{v^2}{2 v_{th}^2})
```
"""
function sample!( pg::ParticleGroup{1,1}, α::Float64, k::Float64, σ::Float64, mesh::OneDGrid)


    s = Sobol.SobolSeq(2)

    nbpart = pg.n_particles

    for i = 1:nbpart
        v = σ * sqrt(-2 * log((i - 0.5) / nbpart))
        r1, r2 = Sobol.next!(s)
        θ = r1 * 2π
        pg.array[1,i] = newton(r2, α, k) * mesh.dimx * k / 2π
        pg.array[2,i] = v * cos(θ) 
        pg.array[3,i] = mesh.dimx
    end

end

