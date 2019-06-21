export ParticleGroup

abstract type AbstractParticleGroup end

mutable struct ParticleGroup{D,V} <:  AbstractParticleGroup 

    dims              :: Tuple{Int64, Int64}
    n_particles       :: Int64
    n_total_particles :: Int64
    particle_array    :: Array{Float64, 2} 
    common_weight     :: Float64
    charge            :: Float64
    mass              :: Float64
    n_weights         :: Int64
    q_over_m          :: Float64

    function ParticleGroup{D,V}( n_particles, 
                                 n_total_particles, 
                                 charge, 
                                 mass, 
                                 n_weights) where {D, V}

        dims = (D, V)
        particle_array = zeros( Float64, (sum(dims)+n_weights, n_particles)) 
        common_weight  = 1.0
        q_over_m = charge / mass

        new( dims,
             n_particles,
             n_total_particles,
             particle_array,
             common_weight,
             charge,
             mass,
             n_weights,
             q_over_m )
    end
end 

export get_x

"""  
Get position
"""
@generated function get_x( p :: ParticleGroup{D,V}, i ) where {D, V}

    :(p.particle_array[1:$D, i])
    
end 

export get_v

"""  
Get velocities
"""
@generated function get_v( p :: ParticleGroup{D,V}, i ) where {D, V}

    :(p.particle_array[$D+1:$D+$V, i])
end


"""
Get charge of particle (q * particle_weight)
"""
@generated function get_charge( p :: ParticleGroup{D,V}, i; i_wi=1) where {D, V}

    :(p.charge * p.particle_array[$D+$V+i_wi, i] * p.common_weight)

end 


"""
Get mass of particle (m * particle_weight)
"""
@generated function get_mass( p :: ParticleGroup{D,V}, i; i_wi=1) where {D,V}

    :(p.mass * p.particle_array[$D+$V+i_wi, i] * p.common_weight)

end

"""
Get particle weights
"""
@generated function get_weights( p :: ParticleGroup{D,V}, i) where {D, V}

    :(p.particle_array[$D+$V+1:$D+$V+p.n_weights, i])

end

"""
Set position of particle @ i
"""
@generated function set_x( p :: ParticleGroup{D,V}, i, x ) where {D, V}

    :(p.particle_array[1:$D, i] .= x)
    
end

"""
Set velocity of particle @ i
"""
@generated function set_v( p :: ParticleGroup{D,V}, i, x ) where {D, V}

    :(p.particle_array[$D+1:$D+$V, i] .= x)
    
end
  
"""
Set weights of particle @ i
"""
@generated function set_weights( p :: ParticleGroup{D,V}, i, x ) where {D, V}

    :(p.particle_array[$D+$V+1:$D+$V+p.n_weights, i] .= x)
    
end

"""
Set the common weight
"""
function set_common_weight( p :: AbstractParticleGroup, x ) 

    p.common_weight = x
    
end
