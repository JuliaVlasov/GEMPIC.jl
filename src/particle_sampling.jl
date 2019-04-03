using Sobol, Random
using VlasovBase
using Distributions

export ParticleSampler

"""
Particle initializer class with various functions to initialize a particle.
"""
struct ParticleSampler

    sampling_type :: Symbol
    dims          :: Tuple{Int64, Int64}
    n_particles   :: Int
    symmetric     :: Bool
    seed          :: Int64

    function ParticleSampler( sampling_type :: Symbol, 
                              symmetric     :: Bool, 
                              dims          :: Tuple{Int64,Int64}, 
                              n_particles   :: Int64)

        if !(sampling_type in [:random, :sobol])
            throw(ArgumentError("Sampling type $sampling_type 
                      not implemented"))
        end

        # Make sure that the particle number is conforming 
        # with symmetric sampling (if necessary)
        if symmetric
            ncopies = 2^(sum(dims))
            np = mod(n_particles,ncopies)
            if np != 0
                n_particles += np
            end
        else
            ncopies = 1
        end

        seed = 1234

        new( sampling_type, dims, n_particles, symmetric, seed)

    end

end
    
export sample

"""
Sample from distribution defined by \a params
- `xmin` : lower bound of the domain
- `Lx`   : length of the domain.
"""
function sample( ps   :: ParticleSampler, 
                 pg   :: ParticleGroup{1,2}, 
                 df   :: CosGaussian, 
                 mesh :: Mesh )

    if ps.symmetric 
        sample_sym( ps, pg, df, mesh )
    else
        sample_all( ps, pg, df, mesh )
    end 
    
end


"""
Helper function for pure sampling
"""
function sample_all( ps, pg, df, mesh )

    ndx, ndv = df.dims

    x = zeros( ndx )
    v = zeros( ndv )

    n_rnds = 0
    if df.n_gaussians > 1
       n_rnds = 1
    end

    delta = zeros(df.n_gaussians)
    for i_v=1:df.n_gaussians
       delta[i_v] = sum(df.delta[1:i_v])
    end
    
    n_rnds += ndx + ndv
    rdn = zeros(ndx + ndv + 1)
    
    # 1/Np in common weight
    set_common_weight(pg, (1.0/pg.n_total_particles))

    if ps.sampling_type == :sobol
       rng_sobol  = SobolSeq(ndx)
    else
       rng_random = MersenneTwister(ps.seed)
    end 

    d = Normal()
   
    for i_part = 1:pg.n_particles

       if ps.sampling_type == :sobol
           x .= mesh.xmin .+ next!(rng_sobol) .* mesh.Lx
       else
           x .= mesh.xmin .+ rand(rng_random, ndx) .* mesh.Lx
       end

       # Set weight according to value of perturbation
       wi = eval_x_density(df, x) .* prod(mesh.Lx)

       v .= rand!(d, v)

       # For multiple Gaussian, draw which one to take
       rnd_no = rdn[ndx+ndv+1]
       i_gauss = 1
       while( rnd_no > delta[i_gauss] )
          i_gauss += 1
       end
       v .= v .* df.v_thermal[:,i_gauss] .+ df.v_mean[:,i_gauss]
       
       # Copy the generated numbers to the particle
       set_x(pg, i_part, x)
       set_v(pg, i_part, v)
       # Set weights.
       set_weights(pg, i_part, wi)

    end
       
    
end

"""
Helper function for antithetic sampling in 1d2v
"""
function sample_sym( ps, pg, df, mesh )

#=
    sll_int32 :: n_rnds
    sll_real64                                         :: x(3),v(3)
    sll_int32                                          :: i_part
    sll_int32                                          :: i_v
    sll_real64, allocatable                            :: rdn(:)
    sll_real64                                         :: wi(1)
    sll_real64                                         :: rnd_no
    sll_int32                                          :: ip, i_gauss
    sll_real64                                         :: delta(params%n_gaussians)
=#

    ndx, ndv = df.dims

    x = zeros(Float64, ndx)
    v = zeros(Float64, ndv)
    
    n_rnds = 0
    (df.n_gaussians > 1) && (n_rnds = 1)

    delta = zeros(Float64, df.n_gaussians)
    for i_v = 1:df.n_gaussians
       delta[i_v] = sum(df.delta[1:i_v])
    end
    
    n_rnds = n_rnds + ndx + ndv
    rdn    = zeros( ndx + ndv + 1 )
    
    # 1/Np in common weight
    set_common_weight(pg, 1.0/pg.n_total_particles)

    if ps.sampling_type == :sobol
       rng_sobol  = SobolSeq(ndx + ndv + 1)
    else
       rng_random = MersenneTwister(ps.seed)
    end 

    dnormal = Normal()

    i_gauss = 1   :: Int64
    wi      = 0.0 :: Float64
    
    for i_part = 1:pg.n_particles

        ip = i_part % 8

        if ip == 1

            # Generate Random or Sobol numbers on [0,1]
            if ps.sampling_type == :sobol
               rdn .= next!(rng_sobol)
            else
               rdn .= rand!(rng_random, rdn)
            end 
            
            # Transform rdn to the interval
            x[1:ndx] .= mesh.xmin .+ mesh.Lx .* rdn[1:ndx]
            
            # Set weight according to value of perturbation
            wi = eval_x_density(df, x[1:ndx] * prod(mesh.Lx))
            
            # Maxwellian distribution of the temperature

            v .= rand!(dnormal, v)

            # For multiple Gaussian, draw which one to take
            rnd_no = rdn[ndx+ndv+1]
            i_gauss = 1
            while i_gauss < df.n_gaussians && rnd_no > delta[i_gauss]
               i_gauss += 1
            end

            v .= v .* df.v_thermal[:,i_gauss] .+ df.v_mean[:,i_gauss]

        elseif ip == 5

            x[1] = mesh.Lx[1] - x[1] + 2.0 * mesh.xmin[1]

        elseif ip % 2 == 0

            v[1] = - v[1] + 2.0 * df.v_mean[1,i_gauss]

        else          

            v[2] = - v[2] + 2.0 * df.v_mean[2,i_gauss]

        end
             
        # Copy the generated numbers to the particle

        set_x( pg, i_part, x)
        set_v( pg, i_part, v)
        set_weights( pg, i_part, wi)
       
    end

end
