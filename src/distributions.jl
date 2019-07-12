export CosGaussian

"""
    CosGaussian( dims, n_cos, n_gaussians, kx, alpha, v_thermal, v_mean, δ )

# Data type for parameters of initial distribution

## Descriptors for various distributions

- sumcos_onegaussian : Descriptor for 
```math
(1+\\sum \\cos( kx * x_i)) * \\exp (-0.5(v-v_mean)^2/v_thermal^2)
```
- cossum_onegaussian : Descriptor for 
```math
(1+ \\cos( \\sum kx_i * x_i)) * \\exp (-0.5(v-v_mean)^2/v_thermal^2)
```
- sumcos_twogaussian : sumcos_onegaussian but with sum of two Gaussians
- cossum_twogaussian : as sumcos_onegaussian but with sum of two Gaussians

## Parameters

- `kx`          : values of the wave numbers (first index dimension, 
                  second index for multiple cosines)
- `alpha`       : strength of perturbations
- `v_thermal`   : variance of the Gaussian ( first index velocity dimension, 
                  second index multiple Gaussians)
- `v_mean`      : mean value of the Gaussian ( first index velocity dimension,
                  second index multiple Gaussians)
- `delta`       : Portion of each Gaussian
- `normal`      : Normalization constant of each Gaussian
- `n_gaussians` : Number of Gaussians
- `n_cos`       : Number of cosines
"""
struct CosGaussian

    dims        :: Tuple{Int64, Int64}
    n_cos       :: Int64
    n_gaussians :: Int64
    kx          :: Array{Float64, 2} 
    alpha       :: Vector{Float64}
    v_thermal   :: Array{Float64, 2}
    v_mean      :: Array{Float64, 2}
    delta       :: Vector{Float64}
    normal      :: Vector{Float64}

    function CosGaussian( dims        :: Tuple{Int64, Int64},
                          n_cos       :: Int64, 
                          n_gaussians :: Int64,
                          kx          :: Array{Float64, 2}, 
                          alpha       :: Vector{Float64}, 
                          v_thermal   :: Array{Float64, 2}, 
                          v_mean      :: Array{Float64, 2}, 
                          δ           :: Float64 )

    

        delta    = zeros(n_gaussians)
        delta[1] = δ
        if n_gaussians > 1
            delta[2] = 1.0 - δ
        end

        normal = zeros(n_gaussians)
        for j in 1:n_gaussians
            normal[j] = 1.0/((2π)^(0.5*dims[2])*prod(v_thermal[:,j]))
        end

        @show dims        
        @show n_cos       
        @show n_gaussians 
        @show kx           
        @show alpha       
        @show v_thermal   
        @show v_mean      
        @show delta       
        @show normal      

        new( dims, n_cos, n_gaussians, kx, alpha, v_thermal, 
             v_mean, delta, normal )


    end

end

export eval_x_density
function eval_x_density( self :: CosGaussian, x )
    
    fval = 1.0
    for j=1:self.n_cos
       fval += self.alpha[j] * cos( sum(self.kx[:,j] .* x) )
    end
    fval

end
  
export eval_v_density
function eval_v_density( self :: CosGaussian, v ) 

    fval = 0.0
    for j=1:self.n_gaussians
       fval += self.normal[j] * self.delta[j] * exp( - 0.5 * 
               sum( ((v .- self.v_mean[:,j]) ./ self.v_thermal[:,j]).^2))
    end
    fval

end 

export eval_xv_density
function eval_xv_density( self :: CosGaussian, x, v ) 

    eval_x_density(self, x) * eval_v_density(self, v)

end 

function (self :: CosGaussian)( x :: Float64, v :: Float64 ) 

    eval_x_density( self, x) * eval_v_density( self, v)

end

function (self :: CosGaussian)( x::Vector{Float64}, v::Vector{Float64} ) 

    f = zeros(Float64, (length(x), length(v)))

    for j in eachindex(v), i in eachindex(x)
       f[i,j] = eval_x_density( self, x[i]) * eval_v_density( self, v[j])
    end

    f
end
