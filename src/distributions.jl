abstract type AbstractCosGaussian end

struct CosGaussianParams

    dims        :: Tuple{Int64, Int64}
    n_cos       :: Int64
    n_gaussians :: Int64
    kx          :: Array{Vector{Float64}, 1}
    alpha       :: Vector{Float64}
    v_thermal   :: Array{Vector{Float64}, 1}
    v_mean      :: Array{Vector{Float64}, 1}
    normal      :: Vector{Float64}
    delta       :: Vector{Float64}

    function CosGaussianParams( dims      :: Tuple{Int64, Int64},
                                kx        :: Array{Vector{Float64}, 1}, 
                                alpha     :: Vector{Float64}, 
                                v_thermal :: Array{Vector{Float64}, 1}, 
                                v_mean    :: Array{Vector{Float64}, 1},
                                delta     :: Vector{Float64} = [1.0]) 

        n_cos = length(kx) 
        @assert n_cos == length(alpha)
        for i in 1:n_cos
            @assert length(kx[i]) == dims[1]
        end
        n_gaussians = length(v_thermal)
        @assert n_gaussians == length(v_mean)
        @assert all( v_thermal .!= 0.0)
        normal = zeros(n_gaussians)
        for j in 1:n_gaussians
            @assert length(v_thermal[j]) == dims[2]
            normal[j] = 1.0/((2π)^(0.5*dims[2])*prod(v_thermal[j]))
        end
        @assert sum(delta) == 1.0

        new( dims, n_cos, n_gaussians, kx, alpha, v_thermal, v_mean, normal, delta )
    end

end

export CosSumGaussian

"""
    CosSumGaussian{D,V}( n_cos, n_gaussians, kx, alpha, v_thermal, v_mean, δ )

Data type for parameters of initial distribution

```math
(1+ \\cos( \\sum_{n_{cos}} kx_i x_i)) 
* \\exp \\big( -\\frac{1}{2} \\sum_{n_{gaussians}} \\frac{(v-v_{j,mean})^2}{v_{j,thermal}^2} \\big)
```

## Parameters

- `kx`          : values of the wave numbers (first index dimension, second index for multiple cosines)
- `alpha`       : strength of perturbations
- `v_thermal`   : variance of the Gaussian (array of  velocity vectors).
- `v_mean`      : mean value of the Gaussian (array of velocity vector).
- `normal`      : Normalization constant of each Gaussian
- `n_gaussians` : Number of Gaussians
- `n_cos`       : Number of cosines

# Example

```math
f(x,i\\mathbf{v},t=0)=\\frac{1}{2\\pi\\sigma_1\\sigma_2} \\exp \\Big( - \\frac{1}{2} \\big( \\frac{v_1^2}{\\sigma_1^2} + \\frac{v_2^2}{\\sigma_2^2} \\big) \\Big) ( 1 + \\alpha \\cos(kx)), \\qquad x \\in [0, 2\\pi / k),
```

```julia
alpha     = [0.01]
kx        = [[1.25]]
v_thermal = [[σ₁],  [σ₂]]
v_mean    = [[0.0], [0.0]]
```


"""
struct CosSumGaussian{D, V} <: AbstractCosGaussian

    dims        :: Tuple{Int64, Int64}
    params      :: CosGaussianParams

    function CosSumGaussian{D, V}( kx        :: Array{Vector{Float64}, 1}, 
                                   alpha     :: Vector{Float64}, 
                                   v_thermal :: Array{Vector{Float64}, 1}, 
                                   v_mean    :: Array{Vector{Float64}, 1},
                                   delta     :: Vector{Float64} = [1.0] ) where {D,V}

        dims   = (D, V)
        params = CosGaussianParams( dims, kx, alpha, v_thermal, v_mean, delta )

        new( dims, params )
    end

end

function eval_x_density( f :: CosSumGaussian, x :: Union{Float64,Vector{Float64}} )
    
    fval = 1.0
    for j=1:f.params.n_cos
       fval += f.params.alpha[j] * cos( sum(f.params.kx[j] .* x) )
    end
    fval

end
  


export SumCosGaussian

"""
    SumCosGaussian( dims, n_cos, n_gaussians, kx, alpha, v_thermal, v_mean, δ )
Data type for parameters of initial distribution
```math
(1+ \\sum_{n_{cos}} \\cos(  kx_i * x_i)) \\exp 
\\big( -\\frac{1}{2} \\sum_{n_{gaussians}} \\frac{(v-v_{mean})^2}{v_{thermal}^2} \\big)
```
## Parameters
- `kx`          : values of the wave numbers (first index dimension, second index for multiple cosines)
- `alpha`       : strength of perturbations
- `v_thermal`   : variance of the Gaussian ( first index velocity dimension, second index multiple Gaussians)
- `v_mean`      : mean value of the Gaussian ( first index velocity dimension, second index multiple Gaussians)
- `normal`      : Normalization constant of each Gaussian
- `n_gaussians` : Number of Gaussians
- `n_cos`       : Number of cosines

# Example

```math
f(x,i\\mathbf{v},t=0)=\\frac{1}{2\\pi\\sigma_1\\sigma_2} \\exp \\Big( - \\frac{1}{2} \\big( \\frac{v_1^2}{\\sigma_1^2} + \\frac{v_2^2}{\\sigma_2^2} \\big) \\Big) ( 1 + \\alpha_1 \\cos(k_1 x) + \\alpha_2 \\cos(k_2 x) ), \\qquad x \\in [0, 2\\pi / k),
```
```julia
kx        = [[k₁],[k₂]]
v_thermal = [[σ₁],[σ₂]]
v_mean    = [[0.0],[0.0]]
alpha     = [α₁, α₂]

df = SumCosGaussian{1,2}( kx, alpha, v_thermal, v_mean )

```
"""
struct SumCosGaussian{D,V}
value :: Int64
    double :: Int64
    dims   :: Tuple{Int64,Int64}
    params :: CosGaussianParams

    function SumCosGaussian{D,V}( kx        :: Array{Array{Float64,1},1}, 
                                  alpha     :: Vector{Float64}, 
                                  v_thermal :: Array{Array{Float64,1},1}, 
                                  v_mean    :: Array{Array{Float64,1},1},
                                  delta     :: Array{Float64,1} = [1.0]
                                  ) where {D,V}

        dims   = (D, V)
        params = CosGaussianParams( dims, kx, alpha, v_thermal, v_mean, delta)
        new( dims, params )
    end

end

function eval_x_density( f :: SumCosGaussian, x :: Union{Float64,Vector{Float64}} )
    
    fval = 1.0
    for j=1:f.params.n_cos
       fval += f.params.alpha[j] * cos( sum(f.params.kx[:,j] .* x) )
    end
    fval

end

function eval_v_density( f :: AbstractCosGaussian, v :: Union{Float64,Vector{Float64}} ) 

    fval = 0.0
    for j=1:f.params.n_gaussians
       fval += f.params.normal[j] * f.params.delta[j] .* exp( - 0.5 * 
               sum( ((v .- f.params.v_mean[:,j]) ./ f.params.v_thermal[j]).^2))
    end
    fval

end 

function( f :: CosSumGaussian )( x, v )

    eval_x_density( f, x) * eval_v_density( f, v)

end

function( f :: SumCosGaussian )( x, v )

    eval_x_density( f, x) * eval_v_density( f, v)

end
