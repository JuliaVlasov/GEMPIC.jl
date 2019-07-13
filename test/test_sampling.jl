@testset "Sampling" begin

function test_sampling( sampling_type :: Symbol, 
                        symmetric     :: Bool, 
                        pg            :: ParticleGroup{D,V}, 
                        df            :: GEMPIC.AbstractCosGaussian ) where {D, V}

   mean  = zeros(3)
   sigma = zeros(3)
   xi    = zeros(3)

   n_particles = pg.n_particles
   
   sampling = ParticleSampler{D,V}( sampling_type, symmetric, n_particles )

   sample!( pg, sampling, df, mesh )
   
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

n_particles = 100000
xmin        = 1.0 :: Float64
xmax        = 4π + 1.0
Lx          = xmax - xmin  
nx          = 64

mesh = Mesh( xmax, xmin, nx)

pg = ParticleGroup{1,2}(n_particles, 1.0, 1.0, 1)

params = ( k = [[0.5]],
           α = [0.01],
           σ = [[0.1, 2.0]],
           μ = [[0.0, 0.0]]
)

df1 = CosSumGaussian{1,2}(params...)

mean_ref  = [Lx*0.5+xmin, 0.0, 0.0]
sigma_ref = [Lx^2/12.0,   params.σ[1][1]^2, params.σ[1][2]^2 ]

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
    k = [[0.5]],
    α = [0.01],
    σ = [[0.1, 2.0], [2.0, 2.0]],
    μ = [[0.0, 0.0], [1.0, 1.0]],
    δ = [0.7, 0.3]
)

df2 = CosSumGaussian{1,2}(params...)

mean_ref = [Lx*0.5+xmin, 0.3, 0.3]
sigma_ref[1] = Lx^2/12.0
for j=1:2
    sigma_ref[j+1] = - (df2.params.δ[1] * df2.params.μ[j][1]
                     + (df2.params.δ[2] 
                     -  df2.params.δ[1]) * df2.params.μ[j][2])^2
    for k=1:2
        sigma_ref[j+1] += df2.params.δ[k] * (df2.params.σ[j][k]^2 + df2.params.μ[j][k]^2)
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
