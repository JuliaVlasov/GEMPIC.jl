export ParticleGroup1D2V

mutable struct ParticleGroup1D2V

    n_particles       :: Int64
    n_total_particles :: Int64
    particle_array    :: Array{Float64, 2} 
    common_weight     :: Float64
    charge            :: Float64
    mass              :: Float64
    n_weights         :: Int64

    function ParticleGroup1D2V( n_particles, 
                            n_total_particles, 
                            charge, 
                            mass, 
                            n_weights)

        particle_array = zeros( Float64, (3+n_weights, n_particles)) 
        common_weight  = 1.0

        new( n_particles,
             n_total_particles,
             particle_array,
             common_weight,
             charge,
             mass,
             n_weights)
    end
end 


"""  
Get position
"""
function get_x( self :: ParticleGroup1D2V, i )

    self.particle_array[1, i]
    
end 

"""  
Get velocities
"""
function get_v( self :: ParticleGroup1D2V, i )

    self.particle_array[2:3, i]
    
end


"""
Get charge of particle (q * particle_weight)
"""
function get_charge( self :: ParticleGroup1D2V, i; i_wi=1)

    self.charge * self.particle_array[3+i_wi, i] * self.common_weight

end 


"""
Get mass of particle (m * particle_weight)
"""
function get_mass( self :: ParticleGroup1D2V, i; i_wi=1) 

    self.mass * self.particle_array[3+i_wi, i] * self.common_weight

end

"""
Get particle weights
"""
function get_weights( self :: ParticleGroup1D2V, i) 

    self.particle_array[4:3+self.n_weights, i]

end

"""
Set position of particle @ i
"""
function set_x( self :: ParticleGroup1D2V, i, x )

    self.particle_array[1, i] = x
    
end

"""
Set velocity of particle @ i
"""
function set_v( self :: ParticleGroup1D2V, i, x )

    self.particle_array[2:3, i] .= x
    
end
  
"""
Set weights of particle @ i
"""
function set_weights( self :: ParticleGroup1D2V, i, x )

    self.particle_array[4:3+self.n_weights, i] .= x
    
end

"""
Set the common weight
"""
function set_common_weight( self :: ParticleGroup1D2V, x )

    self.common_weight = x
    
end
