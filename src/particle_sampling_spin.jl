"""
    sample!( pg, ps, df, mesh)

Sample from a Particle sampler

- `pg`   : Particle group
- `ps`   : Particle sampler
- `df`   : Distribution function
- `xmin` : lower bound of the domain
- `Lx`   : length of the domain.
"""
function sample!( pg   :: ParticleGroup{1,1}, 
                  ps   :: ParticleSampler, 
                  df   :: CosSumGaussianSpin, 
                  mesh :: AbstractGrid )

    s = zeros( pg.n_spin )

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
    
    n_rnds += 2
    rdn = zeros(3)
    
    if ps.sampling_type == :sobol
       rng_sobol  = SobolSeq(1)
    end 

    rng_random = MersenneTwister(ps.seed)

    d = Normal()

    for i_part = 1:(pg.n_particles)  
       
       if ps.sampling_type == :sobol
           x = mesh.xmin + Sobol.next!(rng_sobol)[1] * mesh.dimx
       else
           x = mesh.xmin + rand(rng_random) * mesh.dimx
       end

       v = rand(rng_random, d)

       # For multiple Gaussian, draw which one to take
       rnd_no = rdn[3]
        
       i_gauss = 1
       while( rnd_no > δ[i_gauss] )
          i_gauss += 1
       end
       v = v * df.params.σ[i_gauss] + df.params.μ[i_gauss]

#for a peaked initial condition in the s direction 
       s = [0, 0, 1]
#for a uniformly distributed initial condition on the sphere
#       for tt = 1:10
#            s[1] = randn()
#            s[2] = randn()
#            s[3] = randn()
#            if norm(s) > 10^(-4)
#                break
#            end
#        end
#        s .= s./norm(s)

       # Set weight according to value of perturbation
       w  = eval_x_density(df, x) * mesh.dimx
        
       # Copy the generated numbers to the particle
       set_x!(pg, i_part, x[1])
       set_v!(pg, i_part, v)
       set_spin!(pg, i_part, 1, s[1])
       set_spin!(pg, i_part, 2, s[2])
       set_spin!(pg, i_part, 3, s[3])
       # Set weights.
       set_weights!(pg, i_part, w)
        
    end
       
end
