@testset "Sampling" begin

using VlasovBase

xmin     = 1.0 :: Float64
Lx       = 4π  

function test_sampling( sampling_type :: Symbol, 
                        symmetric     :: Bool, 
                        pg            :: ParticleGroup{D,V}, 
                        df            :: CosGaussian ) where {D, V}

   mean  = zeros(3)
   sigma = zeros(3)
   xi    = zeros(3)

   n_particles = pg.n_particles
   
   sampling = ParticleSampler( sampling_type, symmetric, (D, V), n_particles )

   sample( sampling, pg, df, xmin, Lx )
   
   for i_part = 1:n_particles
       xi = get_x(pg, i_part)
       vi = get_v(pg, i_part)
       mean[1] += xi[1]
       mean[2] += vi[1]
       mean[3] += vi[2]
   end

   mean = mean/n_particles

   for i_part = 1:n_particles
       xi = get_x(pg, i_part)
       vi = get_v(pg, i_part)
       sigma[1] += (xi[1] - mean[1])^2
       sigma[2] += (vi[1] - mean[2])^2
       sigma[3] += (vi[2] - mean[3])^2
   end
   
   sigma = sigma/(n_particles-1)
   
   mean, sigma
   
end

n_particles = 80000

pg = ParticleGroup{1,2}(n_particles, n_particles, 1.0, 1.0, 1)

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

@info "Sobol non symmetric"
mean, sigma = test_sampling( :sobol, false, pg, df1 )
@show abs.(mean  .- mean_ref)
@show abs.(sigma .- sigma_ref)
@test maximum(abs.(mean .- mean_ref)) ≈ 0.0 atol = 1e2/sqrt(n_particles)

@info "Sobol symmetric"
mean, sigma = test_sampling( :sobol, true,  pg, df1 )
@show abs.(mean  .- mean_ref)
@show abs.(sigma .- sigma_ref)
@test maximum(abs.(mean .- mean_ref)) ≈ 0.0 atol = 1e-12

@info "Random non symmetric"
mean, sigma = test_sampling( :random, false, pg, df1)
@show abs.(mean  .- mean_ref)
@show abs.(sigma .- sigma_ref)
@test maximum(abs.(mean .- mean_ref)) ≈ 0.0 atol = 1e2/sqrt(n_particles)

@info "Random symmetric"
mean, sigma = test_sampling( :random, true,  pg, df1)
@show abs.(mean  .- mean_ref)
@show abs.(sigma .- sigma_ref)
@test maximum(abs.(mean .- mean_ref)) ≈ 0.0 atol = 1e-12
   
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
    δ           = 0.7
)

df2 = CosGaussian(params...)

mean_ref = [Lx*0.5+xmin, 0.3, 0.3]
sigma_ref[1] = Lx^2/12.0
for j=1:2
    sigma_ref[j+1] = - (df2.delta[1] * df2.v_mean[j,1]
                      +(df2.delta[2] - df2.delta[1])*df2.v_mean[j,2])^2
    for k=1:2
        sigma_ref[j+1] += df2.delta[k] * (df2.v_thermal[j,k]^2
                                       +  df2.v_mean[j,k]^2)
    end
end
  
@info "Sobol non symmetric"
mean, sigma =  test_sampling( :sobol, false, pg, df2) 
@show abs.(mean  .- mean_ref)
@show abs.(sigma .- sigma_ref)
@test maximum(abs.(mean .- mean_ref)) ≈ 0.0 atol = 1e2/sqrt(n_particles)


@info "Sobol symmetric"
mean, sigma =  test_sampling( :sobol, true,  pg, df2) 
@show abs.(mean  .- mean_ref)
@show abs.(sigma .- sigma_ref)
@test maximum(abs.(mean .- mean_ref)) ≈ 0.0 atol = 1e2/sqrt(n_particles)

end
