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
    sample!( pg, ps, df, mesh)

Sample from a Particle sampler

- `pg`   : Particle group
- `ps`   : Particle sampler
- `df`   : Distribution function
- `xmin` : lower bound of the domain
- `dimx` : length of the domain.
"""
function sample!(
    pg::ParticleGroup{1,2}, ps::ParticleSampler, df::AbstractCosGaussian, mesh::AbstractGrid
)
    if ps.symmetric
        sample_sym(ps, pg, df, mesh)
    else
        sample_all(ps, pg, df, mesh)
    end
end

"""
    sample_all( ps, pg, df, mesh )

Helper function for pure sampling
"""
function sample_all(ps, pg::ParticleGroup{1,2}, df::AbstractCosGaussian, mesh)
    ndx, ndv = df.dims

    x = zeros(ndx)
    v = zeros(ndv)

    n_rnds = 0
    if df.params.n_gaussians > 1
        n_rnds = 1
    end

    δ = zeros(df.params.n_gaussians)
    for i_v in 1:(df.params.n_gaussians)
        δ[i_v] = sum(df.params.δ[1:i_v])
    end

    n_rnds += ndx + ndv
    rdn = zeros(ndx + ndv + 1)

    if ps.sampling_type == :sobol
        rng_sobol = SobolSeq(ndx)
    else
        rng_random = MersenneTwister(ps.seed)
    end

    d = Normal()

    for i_part in 1:(pg.n_particles)
        if ps.sampling_type == :sobol
            x .= mesh.xmin .+ Sobol.next!(rng_sobol) * mesh.dimx
        else
            x .= mesh.xmin .+ rand(rng_random, ndx) * mesh.dimx
        end

        # Set weight according to value of perturbation
        w = eval_x_density(df, x) .* mesh.dimx

        v .= rand!(d, v)

        # For multiple Gaussian, draw which one to take
        rnd_no = rdn[ndx + ndv + 1]
        i_gauss = 1
        while (rnd_no > δ[i_gauss])
            i_gauss += 1
        end
        v .= v .* df.params.σ[i_gauss] .+ df.params.μ[i_gauss]

        # Copy the generated numbers to the particle
        set_x!(pg, i_part, x)
        set_v!(pg, i_part, v)

        # Set weights.
        set_weights!(pg, i_part, w)
    end
end

"""
    sample_sym( ps, pg, df, mesh )

Helper function for antithetic sampling in 1d2v
"""
function sample_sym(ps, pg, df, mesh)
    ndx, ndv = df.dims

    x = zeros(Float64, ndx)
    v = zeros(Float64, ndv)

    n_rnds = 0
    (df.params.n_gaussians > 1) && (n_rnds = 1)

    δ = zeros(Float64, df.params.n_gaussians)
    for i_v in 1:(df.params.n_gaussians)
        δ[i_v] = sum(df.params.δ[1:i_v])
    end

    n_rnds = n_rnds + ndx + ndv
    rdn = zeros(ndx + ndv + 1)

    if ps.sampling_type == :sobol
        rng_sobol = SobolSeq(ndx + ndv + 1)
    end

    rng_random = MersenneTwister(ps.seed)

    dnormal = Normal()

    i_gauss = 1::Int
    wi = 0.0::Float64

    for i_part in 1:(pg.n_particles)
        ip = i_part % 8

        if ip == 1

            # Generate Random or Sobol numbers on [0,1]
            if ps.sampling_type == :sobol
                rdn .= Sobol.next!(rng_sobol)
            else
                rdn .= rand!(rng_random, rdn)
            end

            # Transform rdn to the interval
            x[1:ndx] .= mesh.xmin .+ mesh.dimx .* rdn[1:ndx]

            # Set weight according to value of perturbation
            wi = eval_x_density(df, x[1:ndx]) * mesh.dimx

            # Maxwellian distribution of the temperature

            v .= rand!(rng_random, dnormal, v)

            # For multiple Gaussian, draw which one to take
            rnd_no = rdn[ndx + ndv + 1]
            i_gauss = 1
            while i_gauss < df.params.n_gaussians && rnd_no > δ[i_gauss]
                i_gauss += 1
            end

            v .= v .* df.params.σ[i_gauss] .+ df.params.μ[i_gauss]

        elseif ip == 5
            x[1] = mesh.dimx - x[1] + 2.0 * mesh.xmin

        elseif ip % 2 == 0
            v[1] = -v[1] + 2.0 * df.params.μ[i_gauss][1]

        else
            v[2] = -v[2] + 2.0 * df.params.μ[i_gauss][2]
        end

        # Copy the generated numbers to the particle

        set_x!(pg, i_part, x)
        set_v!(pg, i_part, v)
        set_weights!(pg, i_part, wi)
    end
end


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
        pg.array[1,i] = newton(r2, α, k)
        pg.array[2,i] = v * cos(θ) 
        pg.array[3,i] = mesh.dimx
    end

end

"""
    sample!( pg::ParticleGroup{2,1}, α, k, σ, mesh::OneDGrid)

Sampling from a probability distribution to initialize a Landau damping in
1D1V space.

```math
f_0(x,v,t) = \\frac{n_0}{2π v_{th}^2} ( 1 + \\alpha cos(k_x x)) exp( - \\frac{v^2}{2 v_{th}^2})
```
"""
function sample!( pg::ParticleGroup{1,2}, α::Float64, k::Float64, σ::Float64, mesh::OneDGrid)

    s = Sobol.SobolSeq(2)
    @assert mesh.dimx ≈ 2π / k

    nbpart = pg.n_particles

    for i = 1:nbpart
        v = σ * sqrt(-2 * log((i - 0.5) / nbpart))
        r1, r2 = Sobol.next!(s)
        θ = r1 * 2π
        pg.array[1,i] = newton(r2, α, k)
        pg.array[2,i] = v * cos(θ) 
        pg.array[3,i] = v * sin(θ) 
        pg.array[4,i] = mesh.dimx
    end

end
