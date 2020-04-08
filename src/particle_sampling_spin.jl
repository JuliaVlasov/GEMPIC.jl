using Sobol, Random, Distributions

export ParticleSampler

"""
    ParticleSampler{D,V,S}( sampling_type, symmetric, dims, n_particles)

Particle initializer class with various functions to initialize a particle.

- `sampling_type` : `:random` or `:sobol`
- `symmetric` : `true` or `false`
- `n_particles` : number of particles
"""
struct ParticleSampler{D,V,S}

    sampling_type :: Symbol
    dims          :: Tuple{Int64, Int64, Int64}
    n_particles   :: Int
    symmetric     :: Bool
    seed          :: Int64

    function ParticleSampler{D,V,S}( sampling_type :: Symbol, 
                                   symmetric     :: Bool, 
                                   n_particles   :: Int64,
                                   seed          :: Int64 = 1234) where {D,V,S}

        if !(sampling_type in [:random, :sobol])
            throw(ArgumentError("Sampling type $sampling_type 
                      not implemented"))
        end

        dims = (D, V, S)

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


        new( sampling_type, dims, n_particles, symmetric, seed)

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
function sample!( pg   :: ParticleGroup{1,1,3}, 
                  ps   :: ParticleSampler, 
                  df   :: AbstractCosGaussian, 
                  mesh :: Mesh )

    if ps.symmetric 
        sample_sym( ps, pg, df, mesh )
    else
        sample_all( ps, pg, df, mesh )
    end 
    
end


"""
    sample_all( ps, pg, df, mesh )

Helper function for pure sampling
"""
function sample_all( ps, pg, df :: AbstractCosGaussian, mesh )

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

"""
    sample_sym( ps, pg, df, mesh )

Helper function for antithetic sampling in 1d2v
"""
function sample_sym( ps, pg, df, mesh )

    ndx, ndv, nds = df.dims

    x = zeros(Float64, ndx)
    v = zeros(Float64, ndv)
    
    n_rnds = 0
    (df.params.n_gaussians > 1) && (n_rnds = 1)

    δ = zeros(Float64, df.params.n_gaussians)
    for i_v = 1:df.params.n_gaussians
       δ[i_v] = sum(df.params.δ[1:i_v])
    end
    
    n_rnds = n_rnds + ndx + ndv
    rdn    = zeros( ndx + ndv + 1 )
    
    # 1/Np in common weight
    set_common_weight(pg, 1.0/pg.n_particles)

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
               rdn .= Sobol.next!(rng_sobol)
            else
               rdn .= rand!(rng_random, rdn)
            end 
            
            # Transform rdn to the interval
            x[1:ndx] .= mesh.xmin .+ mesh.Lx .* rdn[1:ndx]
            
            # Set weight according to value of perturbation
            wi = eval_x_density(df, x[1:ndx]) * prod(mesh.Lx)
            
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

            x[1] = mesh.Lx[1] - x[1] + 2.0 * mesh.xmin[1]

        elseif ip % 2 == 0

            v[1] = - v[1] + 2.0 * df.params.μ[i_gauss][1]

        else          

            v[2] = - v[2] + 2.0 * df.params.μ[i_gauss][2]

        end
             
        # Copy the generated numbers to the particle

        set_x( pg, i_part, x)
        set_v( pg, i_part, v)
        set_weights( pg, i_part, wi)
       
    end

end
