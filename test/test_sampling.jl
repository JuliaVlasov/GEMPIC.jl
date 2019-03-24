
@testset "Sampling" begin

using VlasovBase

function test_sampling( sampling_type, symmetric, dims, pg, params, tolerance )


   mean  = zeros(3)
   sigma = zeros(3)
   xi    = zeros(3)

   n_particles = pg.n_particles
   
   sampling = ParticleSampler( sampling_type, symmetric, dims, n_particles )
   sample( sampling, pg, params, xmin, Lx )
   
   for i_part = 1:n_particles
       xi = get_x(pg, i_part)
       vi = get_v(pg, i_part)
       mean[1] += xi
       mean[2] += vi[1]
       mean[3] += vi[2]
   end

   mean = mean/n_particles

   for i_part = 1:n_particles
       xi = get_x(pg, i_part)
       vi = get_v(pg, i_part)
       sigma[1] += (xi    - mean[1])^2
       sigma[2] += (vi[1] - mean[2])^2
       sigma[3] += (vi[2] - mean[3])^2
   end
   
   sigma = sigma/real(n_particles-1, f64)
   
   @show mean  .= mean  .- mean_ref
   @show sigma .= sigma .- sigma_ref
   
   maximum(abs.(mean)) < tolerance
   
end

n_particles = 80000

pg = ParticleGroup{1,2}(n_particles, n_particles, 1.0, 1.0, 1)

xmin     = 1.0
Lx       = 4π

params = ( dims        = (1,2),
           n_cos       = 1,
           n_gaussians = 1,
           kx          = hcat([0.5]),
           alpha       = [0.01],
           v_thermal   = hcat([0.1, 2.0]),
           v_mean      = hcat([0.0, 0.0]),
           δ           = 0.0
)

df1 = CosGaussian(params...)

mean_ref  = [Lx*0.5+xmin, 0.0, 0.0]
sigma_ref = [Lx^2/12.0,   df1.v_thermal[1,1]^2, df1.v_thermal[2,1]^2 ]

#@test test_sampling( :sobol,  false, [1,2], pg, params, 1e2/sqrt(n_particles))
#@test test_sampling( :sobol,  true,  [1,2], pg, params, 1e-12)
#@test test_sampling( :random, false  [1,2], pg, params, 1e2/sqrt(n_particles))
#@test test_sampling( :random, true,  [1,2], pg, params, 1e-12)
   
# Expected mean:
# 2π+1, 0, 0
# Expected variance:
# 16/12 π^2 , 0.01, 4

params = (
    dims        = (1,2),
    n_cos       = 1,
    n_gaussians = 2,
    kx          = hcat([0.5]),
    alpha       = [0.01],
    v_thermal   = hcat([0.1, 2.0], [2.0, 2.0]),
    v_mean      = hcat([0.0, 0.0], [1.0, 1.0]),
    δ           = 0.7,
)

df2 = CosGaussian(params...)

mean_ref = [Lx*0.5+xmin, 0.3, 0.3]
sigma_ref[1] = Lx^2/12.0
for j=1:2
    sigma_ref[j+1] = - (df2.delta[1] * df2.v_mean[j,1]
                      +(df2.delta[2] - df2.delta[1])*df2.v_mean[j,2])^2
    for k=1:2
        sigma_ref[j+1] += df2.delta[k] * (df2.v_thermal[j,k]^2
                                            +df2.v_mean[j,k]^2)
    end
end
  
#@test test_sampling( :sobol, false, [1,2], pg, params, 1e2/sqrt(n_particles))
#@test test_sampling( :sobol, true,  [1,2], pg, params, 1e2/sqrt(n_particles))

end
