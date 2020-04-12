using Printf, FileIO, JLD2

export ParticleGroup

abstract type AbstractParticleGroup end

"""
    ParticleGroup{D,V}( n_particles, charge, mass, q, weights)

- `n_particles` : number of particles 
- `charge`      : charge of the particle species
- `mass`        : mass of the particle species
- `n_weights`   : number of differents weights
"""
struct ParticleGroup{D,V} <:  AbstractParticleGroup 

    dims              :: Tuple{Int, Int}
    n_particles       :: Int
    particle_array    :: Array{Float64, 2} 
    common_weight     :: Ref{Float64}
    charge            :: Float64
    mass              :: Float64
    n_weights         :: Int
    q_over_m          :: Float64
    n_spin            :: Int

    function ParticleGroup{D,V}( n_particles; 
                                 charge = 1.0, 
                                 mass = 1.0, 
                                 n_weights = 1, 
                                 n_spin = 0) where {D, V}

        dims = (D, V)
        particle_array = zeros( Float64, (sum(dims)+n_weights+n_spin, n_particles)) 
        common_weight = 1.0
        q_over_m = charge / mass

        new( dims, n_particles, particle_array, common_weight, charge,
             mass, n_weights, q_over_m, n_spin )
    end
end 

"""  
    get_x( p, i )

Get position of ith particle of p
"""
@inline get_x( p :: ParticleGroup{D,V}, i :: Int ) where {D, V} = p.particle_array[1:D, i]

"""  
    get_v( p, i )

Get velocity of ith particle of p
"""
@inline get_v( p :: ParticleGroup{D,V}, i  :: Int) where {D, V} = p.particle_array[D+1:D+V, i]

"""
    get_charge( p, i; i_wi=1)

Get charge of ith particle of p (q * particle_weight)
"""
@inline function get_charge( p :: ParticleGroup{D,V}, i :: Int; i_wi=1) where {D, V} 

    p.charge * p.particle_array[D+V+i_wi, i] * p.common_weight[]

end

"""
    get_mass( p, i; i_wi=1)

Get mass of ith particle of p (m * particle_weight)
"""
@inline function get_mass( p :: ParticleGroup{D,V}, i :: Int; i_wi=1) where {D,V}

	p.mass * p.particle_array[D+V+i_wi, i] * p.common_weight[]

end

"""
    get_weights( p, i)

Get ith particle weights of group p
"""
@inline function get_weights( p :: ParticleGroup{D,V}, i :: Int) where {D, V}

    p.particle_array[D+V+1:D+V+p.n_weights, i]

end

"""
    get_spin( p, i, j)

Get the jth weight of the ith particle weights of group p
"""
@inline function get_spin( p :: ParticleGroup{D,V}, i::Int, j::Int) where {D, V}

    p.particle_array[D+V+p.n_weights+j, i]

end

"""
    set_spin( p, i, j)

Set the jth weight of the ith particle weights of group p
"""
@inline function set_spin( p :: ParticleGroup{D,V}, i::Int, j::Int, s) where {D, V}

    p.particle_array[D+V+p.n_weights+j, i] = s

end

"""
    set_x( p, i, x ) 

Set position of ith particle of p to x 
"""
@inline function set_x( p :: ParticleGroup{D,V}, i :: Int, x :: Vector{Float64} ) where {D, V}

    for j in 1:D p.particle_array[j, i] = x[j] end
    
end

"""
    set_x( p, i, x)

Set position of ith particle of p to x

!!! note
    if `x` is a scalar value, only the first x dimension will be set.
"""
@inline function set_x( p :: ParticleGroup{D,V}, i :: Int, x :: Float64 ) where {D, V}

    p.particle_array[1, i] = x

end
    

"""
    set_v( p, i, v)

Set velocity of ith particle of p to v
"""
@inline function set_v( p :: ParticleGroup{D,V}, i :: Int, v :: Vector{Float64} ) where {D, V}

    for j in 1:V p.particle_array[D+j, i] = v[j] end
    
end

"""
    set_v( p, i, v)

Set velocity of ith particle of p to v
"""
@inline function set_v( p :: ParticleGroup{D,V}, i :: Int, v :: Float64 ) where {D, V}

    p.particle_array[D+1, i] = v
    
end
  
"""
    set_weights( p, i, w) 

Set weights of ith particle of p to w
"""
function set_weights( p :: ParticleGroup{D,V}, i :: Int, w :: Vector{Float64} ) where {D, V}

    for j in 1:p.n_weights p.particle_array[D+V+j, i] = w[j] end
    
end

"""
    set_weights( p, i, w) 

Set weights of particle @ i
"""
function set_weights( p :: ParticleGroup{D,V}, i :: Int, w :: Float64 ) where {D, V}

    p.particle_array[D+V+1, i] = w
    
end

"""
    set_common_weight( p, x ) 

Set the common weight
"""
function set_common_weight( p :: AbstractParticleGroup, x :: Float64 ) 

    p.common_weight[] = x
    
end

function save( file, step, p :: ParticleGroup{D,V}) where {D, V}

    datafile = @sprintf("%s-%06d.jld2", file, step)

    if p.n_spin > 0

        FileIO.save(datafile, Dict("x" => p.particle_array[1:D,:],
                                   "v" => p.particle_array[D+1:D+V,:], 
                                   "w" => p.particle_array[D+V+1:D+V+p.n_weights,:],
                                   "s" => p.particle_array[D+V+p.n_weights+1:end,:]))

    else

        FileIO.save(datafile, Dict("x" => p.particle_array[1:D,:],
                                   "v" => p.particle_array[D+1:D+V,:],
                                   "w" => p.particle_array[D+V+1:end,:]))
    end

end
