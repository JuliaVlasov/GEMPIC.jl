# ---
# jupyter:
#   jupytext:
#     comment_magics: false
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.1
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

# +
mutable struct ParticleGroup{D, V}
    
    n    :: Int
    dims :: Tuple{Int, Int}
    data :: Array{Float64, 2}
    
    function ParticleGroup{D, V}( n ) where {D, V}
        
        dims = (D, V)
        data = zeros(Float64, (D+V,n))
        
        new(n, dims, data)
        
    end

end

# +
@generated function get_x( p ::ParticleGroup{D,V}, i) where {D, V}
    
    return :(p.data[1:$D, i])

end

# +
@generated function get_v( p ::ParticleGroup{D,V}, i) where {D, V}
    
    return :(p.data[$D+1:$D+$V, i])

end
# -

p1d1v = ParticleGroup{1,1}(5)
p1d2v = ParticleGroup{1,2}(5)

get_x(p1d2v, 1), get_v(p1d2v, 1)

get_x(p1d1v, 1), get_v(p1d1v, 1)


