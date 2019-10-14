export ParticleGroup

abstract type AbstractParticleGroup end

"""

    ParticleGroup{D,V}( n_particles, 
                        charge, 
                        mass, 
                        n_weights)

- `n_particles` : number of particles 
- `charge`      : charge of one particle
- `mass`        : mass of one particle
- `n_weights`   : number of differents weights
"""
mutable struct ParticleGroup{D,V} <:  AbstractParticleGroup 

    dims              :: Tuple{Int64, Int64}
    n_particles       :: Int64
    particle_array    :: Array{Float64, 2} 
    common_weight     :: Float64
    charge            :: Float64
    mass              :: Float64
    n_weights         :: Int64
    q_over_m          :: Float64

    function ParticleGroup{D,V}( n_particles, 
                                 charge, 
                                 mass, 
                                 n_weights) where {D, V}

        dims = (D, V)
        particle_array = zeros( Float64, (sum(dims)+n_weights, n_particles)) 
        common_weight  = 1.0
        q_over_m = charge / mass

        new( dims, n_particles, particle_array, common_weight, charge,
             mass, n_weights, q_over_m )
    end
end 

export get_x

"""  
    get_x( p, i )

Get position of ith particle of p
"""
@generated function get_x( p :: ParticleGroup{D,V}, i :: Int64 ) where {D, V}

    :(p.particle_array[1:$D, i])
    
end 

export get_v

"""  
    get_v( p, i )

Get velocity of ith particle of p
"""
@generated function get_v( p :: ParticleGroup{D,V}, i  :: Int64) where {D, V}

    :(p.particle_array[$D+1:$D+$V, i])
end


"""
    get_charge( p, i; i_wi=1)

Get charge of ith particle of p (q * particle_weight)
"""
@generated function get_charge( p :: ParticleGroup{D,V}, i :: Int64; i_wi=1) where {D, V}

    :(p.charge * p.particle_array[$D+$V+i_wi, i] * p.common_weight)

end 


"""
    get_mass( p, i; i_wi=1)

Get mass of ith particle of p (m * particle_weight)
"""
@generated function get_mass( p :: ParticleGroup{D,V}, i :: Int64; i_wi=1) where {D,V}

    :(p.mass * p.particle_array[$D+$V+i_wi, i] * p.common_weight)

end

"""
    get_weights( p, i)

Get ith particle weights of group p
"""
@generated function get_weights( p :: ParticleGroup{D,V}, i :: Int64) where {D, V}

    :(p.particle_array[$D+$V+1:$D+$V+p.n_weights, i])

end

"""
    set_x( p, i, x ) 

Set position of ith particle of p to x 
"""
@generated function set_x( p :: ParticleGroup{D,V}, i :: Int64, x :: Vector{Float64} ) where {D, V}

    :(for j in 1:$D p.particle_array[j, i] = x[j] end)
    
end

"""
    set_x( p, i, x)

Set position of ith particle of p to x

!!! note
    if `x` is a scalar value, only the first x dimension will be set.
"""
@generated function set_x( p :: ParticleGroup{D,V}, i :: Int64, x :: Float64 ) where {D, V}

    :(p.particle_array[1, i] = x)

end
    

"""
    set_v( p, i, v)

Set velocity of ith particle of p to v
"""
@generated function set_v( p :: ParticleGroup{D,V}, i :: Int64, v :: Vector{Float64} ) where {D, V}

    :(for j in 1:$V p.particle_array[$D+j, i] = v[j] end)
    
end

"""
    set_v( p, i, v)

Set velocity of ith particle of p to v
"""
@generated function set_v( p :: ParticleGroup{D,V}, i :: Int64, v :: Float64 ) where {D, V}

    :(p.particle_array[$D+1, i] = v)
    
end
  
"""
    set_weights( p, i, w) 

Set weights of ith particle of p to w
"""
@generated function set_weights( p :: ParticleGroup{D,V}, i :: Int64, w :: Vector{Float64} ) where {D, V}

    :(for j in 1:p.n_weights p.particle_array[$D+$V+j, i] = w[j] end)
    
end

"""
    set_weights( p, i, w) 

Set weights of particle @ i
"""
@generated function set_weights( p :: ParticleGroup{D,V}, i :: Int64, w :: Float64 ) where {D, V}

    :(p.particle_array[$D+$V+1, i] = w)
    
end

"""
    set_common_weight( p, x ) 

Set the common weight
"""
function set_common_weight( p :: AbstractParticleGroup, x :: Float64 ) 

    p.common_weight = x
    
end
