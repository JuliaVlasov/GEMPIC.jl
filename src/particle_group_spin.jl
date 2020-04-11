using JLD2, FileIO, Printf


export SpinParticleGroup

"""
    ParticleGroup{D,V,S}( n_particles, 
                        charge, 
                        mass, 
                        n_weights)

- `n_particles` : number of particles 
- `charge`      : charge of the particle species
- `mass`        : mass of the particle species
- `n_weights`   : number of differents weights
"""
struct SpinParticleGroup{D,V,S} <:  AbstractParticleGroup 

    dims              :: Tuple{Int64, Int64, Int64}
    n_particles       :: Int64
    particle_array    :: Array{Float64, 2} 
    common_weight     :: Ref{Float64}
    charge            :: Float64
    mass              :: Float64
    n_weights         :: Int64
    q_over_m          :: Float64

    function SpinParticleGroup{D,V,S}( n_particles, 
                                 charge, 
                                 mass, 
                                 n_weights) where {D, V, S}

        dims = (D, V, S)
        particle_array = zeros( Float64, (sum(dims)+n_weights, n_particles)) 
        common_weight  = 1.0
        q_over_m = charge / mass

        new( dims, n_particles, particle_array, common_weight, charge,
             mass, n_weights, q_over_m )
    end
end 


"""  
    get_x( p, i )

Get position of ith particle of p
"""
@inline function get_x( p :: SpinParticleGroup{D,V,S}, i :: Int64 ) where {D, V, S}

    p.particle_array[1:D, i]
    
end 


"""  
    get_v( p, i )

Get velocity of ith particle of p
"""
@inline function get_v( p :: SpinParticleGroup{D,V,S}, i  :: Int64) where {D, V, S}

    p.particle_array[D+V, i]
end


"""  
    get_s1( p, i )

Get s1 of ith particle of p
"""
@inline function get_s1( p :: SpinParticleGroup{D,V,S}, i  :: Int64) where {D, V, S}

    p.particle_array[D+V+1, i]
end


"""  
    get_s2( p, i )

Get s2 of ith particle of p
"""
@inline function get_s2( p :: SpinParticleGroup{D,V,S}, i  :: Int64) where {D, V, S}

    p.particle_array[D+V+2, i]
end


"""  
get_s3( p, i )

Get velocity of ith particle of p
"""
@inline function get_s3( p :: SpinParticleGroup{D,V,S}, i  :: Int64) where {D, V, S}

    p.particle_array[D+V+3, i]
end




"""
    get_charge( p, i; i_wi=1)

Get charge of ith particle of p (q * particle_weight)
"""
@inline function get_charge( p :: SpinParticleGroup{D,V,S}, i :: Int64; i_wi=1) where {D, V, S}

    p.charge * p.particle_array[D+V+S+i_wi, i] * p.common_weight[]

end 


"""
    get_mass( p, i; i_wi=1)

Get mass of ith particle of p (m * particle_weight)
"""
@inline function get_mass( p :: SpinParticleGroup{D,V,S}, i :: Int64; i_wi=1) where {D,V,S}

    p.mass * p.particle_array[D+V+S+i_wi, i] * p.common_weight[]

end

"""
    get_weights( p, i)

Get ith particle weights of group p
"""
@inline function get_weights( p :: SpinParticleGroup{D,V,S}, i :: Int64) where {D, V, S}

    p.particle_array[D+V+S+1, i]

end

"""
    set_x( p, i, x ) 

Set position of ith particle of p to x 
"""
@inline function set_x( p :: SpinParticleGroup{D,V,S}, i :: Int64, x :: Vector{Float64} ) where {D, V, S}

    for j in 1:D p.particle_array[j, i] = x[j] end
    
end

"""
    set_x( p, i, x)

Set position of ith particle of p to x

!!! note
    if `x` is a scalar value, only the first x dimension will be set.
"""
@inline function set_x( p :: SpinParticleGroup{D,V,S}, i :: Int64, x :: Float64 ) where {D, V, S}

    p.particle_array[1, i] = x

end
    

"""
    set_v( p, i, v)

Set velocity of ith particle of p to v
"""
@inline function set_v( p :: SpinParticleGroup{D,V,S}, i :: Int64, v :: Vector{Float64} ) where {D, V, S}

    for j in 1:V p.particle_array[D+j, i] = v[j] end
    
end

"""
    set_v( p, i, v)

Set velocity of ith particle of p to v
"""
@inline function set_v( p :: SpinParticleGroup{D,V,S}, i :: Int64, v :: Float64 ) where {D, V, S}

    p.particle_array[D+V, i] = v
    
end

"""
    set_s1( p, i, v)

Set velocity of ith particle of p to v
"""
@inline function set_s1( p :: SpinParticleGroup{D,V,S}, i :: Int64, v :: Float64 ) where {D, V, S}

    p.particle_array[D+V+1, i] = v
    
end


"""
    set_s2( p, i, v)

Set velocity of ith particle of p to v
"""
@inline function set_s2( p :: SpinParticleGroup{D,V,S}, i :: Int64, v :: Float64 ) where {D, V, S}

    p.particle_array[D+V+2, i] = v
    
end

"""
    set_s3( p, i, v)

Set velocity of ith particle of p to v
"""
@inline function set_s3( p :: SpinParticleGroup{D,V,S}, i :: Int64, v :: Float64 ) where {D, V, S}

    p.particle_array[D+V+3, i] = v
    
end

  
"""
    set_weights( p, i, w) 

Set weights of ith particle of p to w
"""
function set_weights( p :: SpinParticleGroup{D,V,S}, i :: Int64, w :: Vector{Float64} ) where {D, V, S}

    for j in 1:p.n_weights p.particle_array[D+V+S+j, i] = w[j] end
    
end

"""
    set_weights( p, i, w) 

Set weights of particle @ i
"""
function set_weights( p :: SpinParticleGroup{D,V,S}, i :: Int64, w :: Float64 ) where {D, V, S}

    p.particle_array[D+V+S+1, i] = w
    
end


function save( file, step, p :: SpinParticleGroup{D,V,S}) where {D, V, S}

    datafile = @sprintf("%s-%06d.jld2", file, step)

    FileIO.save(datafile, Dict("x" => p.particle_array[1:D,:],
                               "v" => p.particle_array[D+1:D+V,:], 
                               "s" => p.particle_array[D+V+1:D+V+S,:],
                               "w" => p.particle_array[D+V+S+1,:]))

end
