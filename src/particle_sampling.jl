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
        sampling_type::Symbol,
        symmetric::Bool,
        n_particles::Int,
        seed::Int = 1234,
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

        new(sampling_type, dims, n_particles, symmetric, seed)

    end

    function ParticleSampler{D,V}(
        sampling_type::Symbol,
        n_particles::Int,
        seed::Int = 1234,
    ) where {D,V}

        if !(sampling_type in [:random, :sobol])
            throw(ArgumentError("Sampling type $sampling_type 
                      not implemented"))
        end

        dims = (D, V)
        symmetric = false

        new(sampling_type, dims, n_particles, symmetric, seed)

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
    pg::ParticleGroup{1,2},
    ps::ParticleSampler,
    df::AbstractCosGaussian,
    mesh::AbstractGrid,
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
    for i_v = 1:df.params.n_gaussians
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

    for i_part = 1:pg.n_particles
        if ps.sampling_type == :sobol
            x .= mesh.xmin .+ Sobol.next!(rng_sobol) * mesh.dimx
        else
            x .= mesh.xmin .+ rand(rng_random, ndx) * mesh.dimx
        end

        # Set weight according to value of perturbation
        w = eval_x_density(df, x) .* mesh.dimx

        v .= rand!(d, v)

        # For multiple Gaussian, draw which one to take
        rnd_no = rdn[ndx+ndv+1]
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
    for i_v = 1:df.params.n_gaussians
        δ[i_v] = sum(df.params.δ[1:i_v])
    end

    n_rnds = n_rnds + ndx + ndv
    rdn = zeros(ndx + ndv + 1)

    if ps.sampling_type == :sobol
        rng_sobol = SobolSeq(ndx + ndv + 1)
    else
        rng_random = MersenneTwister(ps.seed)
    end

    dnormal = Normal()

    i_gauss = 1::Int
    wi = 0.0::Float64

    for i_part = 1:pg.n_particles
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

            v .= rand!(dnormal, v)

            # For multiple Gaussian, draw which one to take
            rnd_no = rdn[ndx+ndv+1]
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
