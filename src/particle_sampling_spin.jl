using Sobol, Random, Distributions

export SpinParticleSampler

"""
    SpinParticleSampler{D,V,S}( sampling_type, symmetric, dims, n_particles)

Particle initializer class with various functions to initialize a particle.

- `sampling_type` : `:random` or `:sobol`
- `symmetric` : `true` or `false`
- `n_particles` : number of particles
"""
struct SpinParticleSampler{D,V,S}

    sampling_type :: Symbol
    dims          :: Tuple{Int64, Int64, Int64}
    n_particles   :: Int
    seed          :: Int64

    function SpinParticleSampler{D,V,S}( sampling_type :: Symbol, 
                                   n_particles   :: Int64,
                                   seed          :: Int64 = 1234) where {D,V,S}

        if !(sampling_type in [:random, :sobol])
            throw(ArgumentError("Sampling type $sampling_type 
                      not implemented"))
        end

        dims = (D, V, S)

        new( sampling_type, dims, n_particles, seed)

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
- `Lx`   : length of the domain.
"""
function sample!( pg   :: SpinParticleGroup{1,1,3}, 
                  ps   :: SpinParticleSampler, 
                  df   :: AbstractCosGaussian, 
                  mesh :: Mesh )

    ndx, ndv, nds = df.dims

    x = zeros( ndx )
    v = zeros( ndv )
    s = zeros( nds )
    theta = 0.0
    phi = 0.0
    n_rnds = 0
    if df.params.n_gaussians > 1
       n_rnds = 1
    end

    δ = zeros(df.params.n_gaussians)
    for i_v=1:df.params.n_gaussians
       δ[i_v] = sum(df.params.δ[1:i_v])
    end
    
    n_rnds += ndx + ndv
    rdn = zeros(ndx + ndv + 1)
    
    # 1/Np in common weight
    set_common_weight(pg, (1.0/pg.n_particles))

    if ps.sampling_type == :sobol
       rng_sobol  = SobolSeq(ndx)
    else
       rng_random = MersenneTwister(ps.seed)
    end 

    d = Normal()
    #nnn = Int64(pg.n_particles/2)
    for i_part = 1:(pg.n_particles)  
     #for i_part = 1:nnn  #  这里devided by 2
       
       if ps.sampling_type == :sobol
           x .= mesh.xmin .+ Sobol.next!(rng_sobol) .* mesh.Lx
       else
           x .= mesh.xmin .+ rand(rng_random, ndx) .* mesh.Lx
       end

       

       v .= rand!(d, v)

       # For multiple Gaussian, draw which one to take
       rnd_no = rdn[ndx+ndv+1]
        
       i_gauss = 1
       while( rnd_no > δ[i_gauss] )
          i_gauss += 1
       end
       v .= v .* df.params.σ[i_gauss] .+ df.params.μ[i_gauss]
       #
	#=        
        for tt = 1:10
            s[1] = randn()
            s[2] = randn()
            s[3] = randn()
            if norm(s) > 10^(-4)
                break
            end
        end
        
        s .= s./norm(s)
        =#
	s[3] = 1
	s[2] = 0
	s[1] = 0
        # Set weight according to value of perturbation
        w  = eval_x_density(df, x) .* prod(mesh.Lx) 
        
       # Copy the generated numbers to the particle
       set_x(pg, i_part, x)
       set_v(pg, i_part, v)
       set_s1(pg, i_part, s[1])
       set_s2(pg, i_part, s[2])
       set_s3(pg, i_part, s[3])
       # Set weights.
       set_weights(pg, i_part, w)
        
        ############################################################# first loop
       #= 
        if ps.sampling_type == :sobol
           x .= mesh.xmin .+ Sobol.next!(rng_sobol) .* mesh.Lx
       else
           x .= mesh.xmin .+ rand(rng_random, ndx) .* mesh.Lx
       end

       # Set weight according to value of perturbation
       w  = eval_x_density(df, x) .* prod(mesh.Lx)

       v .= rand!(d, v)

       # For multiple Gaussian, draw which one to take
       rnd_no = rdn[ndx+ndv+1]
        
       i_gauss = 1
       while( rnd_no > δ[i_gauss] )
          i_gauss += 1
       end
       v .= v .* df.params.σ[i_gauss] .+ df.params.μ[i_gauss]
       #
        for tt = 1:10
            s[1] = randn()
            s[2] = randn()
            s[3] = randn()
            if norm(s) > 10^(-4)
                break
            end
        end
        
        s .= s./norm(s)
       # Copy the generated numbers to the particle
       set_x(pg, 2*i_part-1, x)
       set_v(pg, 2*i_part-1, v)
       set_s1(pg, 2*i_part-1, s[1])
       set_s2(pg, 2*i_part-1, s[2])
       set_s3(pg, 2*i_part-1, s[3])
       # Set weights.
       set_weights(pg, 2*i_part-1, w)
        ############################################################# second loop
        if ps.sampling_type == :sobol
           x .= mesh.xmin .+ Sobol.next!(rng_sobol) .* mesh.Lx
       else
           x .= mesh.xmin .+ rand(rng_random, ndx) .* mesh.Lx
       end

       # Set weight according to value of perturbation
       w  = eval_x_density(df, x) .* prod(mesh.Lx)

       v .= rand!(d, v)

       # For multiple Gaussian, draw which one to take
       rnd_no = rdn[ndx+ndv+1]
        
       i_gauss = 1
       while( rnd_no > δ[i_gauss] )
          i_gauss += 1
       end
       v .= v .* df.params.σ[i_gauss] .- df.params.μ[i_gauss]
       #
        for tt = 1:10
            s[1] = randn()
            s[2] = randn()
            s[3] = randn()
            if norm(s) > 10^(-4)
                break
            end
        end
        
        s .= s./norm(s)
                
       # Copy the generated numbers to the particle
       set_x(pg, 2*i_part, x)
       set_v(pg, 2*i_part, v)
       set_s1(pg, 2*i_part, s[1])
       set_s2(pg, 2*i_part, s[2])
       set_s3(pg, 2*i_part, s[3])
       # Set weights.
       set_weights(pg, 2*i_part, w)
       =#
    end
       
    
    
    
end

