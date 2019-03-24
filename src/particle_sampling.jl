using Sobol, Random

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
function sample( self           :: ParticleSampler, 
                 particle_group :: ParticleGroup{1,2}, 
                 df             :: CosGaussian, 
                 xmin           :: Float64, 
                 Lx             :: Float64 )

    if self.symmetric 
        sample_sym_1d2v( self, particle_group, params, xmin, Lx )
    else
        sample_all( self, particle_group, params, xmin, Lx )
    end 
    
end


"""
Helper function for pure sampling
"""
function sample_particle_sampling_all( self,  particle_group, params, xmin, Lx )

#=
    sll_int32 :: n_rnds
    sll_real64                                         :: x(3),v(3)
    sll_int32                                          :: i_part
    sll_int32                                          :: i_v, i_gauss
    sll_real64, allocatable                            :: rdn(:)
    sll_real64                                         :: wi(particle_group%n_weights)
    sll_real64                                         :: rnd_no
    sll_real64                                         :: delta(params%n_gaussians)
=#

    n_rnds = 0
    if params.n_gaussians > 1
       n_rnds = 1
    end

    for i_v=1:params.n_gaussians
       delta[i_v] = sum(params.delta[1:i_v])
    end
    
    n_rnds += params.dims[1] + params.dims[2]
    rdn = zeros(params.dims[1]+params.dims[2]+1)
    
    # 1/Np in common weight
    set_common_weight(particle_group, (1.0/particle_group.n_total_particles))


    if self.sampling_type == :sobol
       rng_sobol  = SobolSeq(1)
    else
       rng_random = rand(MersenneTwister(self.seed))
    end 
   
#=
   
    for i_part = 1:particle_group.n_particles

       if self.sampling_type == :sobol
           x = next!(rng_sobol)
       else
           x = rand(rng_random)
       end
       x = xmin + Lx * x

       # Set weight according to value of perturbation
       wi = params.eval_x_density(x) * prod(Lx)

       ! Maxwellian distribution of the temperature
       do i_v = 1,params%dims(2)
          call sll_s_normal_cdf_inv( rdn(i_v+params%dims(1)), 0.0_f64, 1.0_f64, &
               v(i_v))
       end do
       ! For multiple Gaussian, draw which one to take
       rnd_no = rdn(params%dims(1)+params%dims(2)+1)
       i_gauss = 1
       do while( rnd_no > delta(i_gauss) )
          i_gauss = i_gauss+1
       end do
       v(1:params%dims(2)) = v(1:params%dims(2)) * params%v_thermal(:,i_gauss) + params%v_mean(:,i_gauss)
       
       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            wi)
       
    end do
       
    
=#

end

"""
Helper function for antithetic sampling in 1d2v
"""
function sample_sym_1d2v( self, particle_group, params, xmin, Lx )

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

    n_rnds = 0
    if ( params%n_gaussians > 1 ) then
       n_rnds = 1
    end if

    do i_v=1,params%n_gaussians
       delta(i_v) = sum(params%delta(1:i_v))
    end do
    
    n_rnds = n_rnds+params%dims(1)+params%dims(2)
    allocate( rdn(params%dims(1)+params%dims(2)+1) )
    rdn = 0.0_f64
    
    ! 1/Np in common weight
    call particle_group%set_common_weight &
         (1.0_f64/real(particle_group%n_total_particles, f64))

    if ( self%random_numbers == sll_p_random_numbers ) then
       call random_seed(put=self%random_seed)
    end if

    do i_part = 1, particle_group%n_particles
       ip = modulo(i_part, 8 )
       if ( ip == 1) then
          ! Generate Random or Sobol numbers on [0,1]
          select case( self%random_numbers )
          case( sll_p_sobol_numbers )
             call sll_s_i8_sobol( int(n_rnds,8), self%sobol_seed, rdn(1:n_rnds))
          case( sll_p_random_numbers )
             call random_number( rdn(1:n_rnds) )
          end select
          
          ! Transform rdn to the interval
          x(1:params%dims(1)) = xmin + Lx * rdn(1:params%dims(1))
          
          ! Set weight according to value of perturbation
          wi(1) = params%eval_x_density(x(1:params%dims(1)))*product(Lx)
          
          ! Maxwellian distribution of the temperature
          do i_v = 1,params%dims(2)
             call sll_s_normal_cdf_inv( rdn(i_v+params%dims(1)), 0.0_f64, 1.0_f64, &
                  v(i_v))
          end do
          ! For multiple Gaussian, draw which one to take
          rnd_no = rdn(params%dims(1)+params%dims(2)+1)
          i_gauss = 1
          do while( rnd_no > delta(i_gauss) )
             i_gauss = i_gauss+1
          end do
          v(1:params%dims(2)) = v(1:params%dims(2)) * params%v_thermal(:,i_gauss) + params%v_mean(:,i_gauss)
       elseif ( ip == 5 ) then
          x(1) = Lx(1) - x(1) + 2.0_f64*xmin(1)
       elseif ( modulo(ip,2) == 0 ) then
          v(1) = -v(1) + 2.0_f64*params%v_mean(1,i_gauss)
       else          
          v(2) = -v(2) + 2.0_f64*params%v_mean(2,i_gauss)
       end if
          
       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            wi)
       
    end do    
=#

end
  

