abstract type AbstractCosGaussian end

"""
    CosGaussianParams( dims, k, α, σ, μ, δ ) 

Parameters of a distribution with is a product of a Cosine 
distribution along x and a Normal distribution along v.

- `n_gaussians` : Number of Gaussians
- `n_cos`       : Number of cosines
- `normal`      : Normalization constant of each Gaussian
"""
struct CosGaussianParams

    dims        :: Tuple{Int64, Int64, Int64}
    n_cos       :: Int64
    n_gaussians :: Int64
    k           :: Array{Vector{Float64}, 1}
    α           :: Vector{Float64}
    σ           :: Array{Vector{Float64}, 1}
    μ           :: Array{Vector{Float64}, 1}
    normal      :: Vector{Float64}
    δ           :: Vector{Float64}

    function CosGaussianParams( dims :: Tuple{Int64, Int64, Int64},
                                k    :: Array{Vector{Float64}, 1}, 
                                α    :: Vector{Float64}, 
                                σ    :: Array{Vector{Float64}, 1}, 
                                μ    :: Array{Vector{Float64}, 1},
                                δ    :: Vector{Float64} = [1.0]) 

        n_cos = length(k) 
        @assert n_cos == length(α)
        for i in 1:n_cos
            @assert length(k[i]) == dims[1]
        end
        n_gaussians = length(σ)
        @assert n_gaussians == length(μ)
        @assert all( σ .!= 0.0)
        normal = zeros(n_gaussians)
        for j in 1:n_gaussians
            @assert length(σ[j]) == dims[2]
            normal[j] = 1.0/((2π)^(0.5*dims[2])*prod(σ[j]))
        end
        @assert sum(δ) == 1.0

        new( dims, n_cos, n_gaussians, k, α, σ, μ, normal, δ )
    end

end

export CosSumGaussian

"""
    CosSumGaussian{D,V}( n_cos, n_gaussians, k, α, σ, μ, δ )

Data type for parameters of initial distribution

```math
(1+ \\cos( \\sum^{n_{cos}}_{i=1} k_i x)) 
\\cdot 
\\sum_{j=1}^{n_{gaussians}} 
\\delta_j 
\\exp \\big( -\\frac{1}{2} 
\\frac{(v-\\mu_j)^2}{\\sigma_j^2} \\big)
```

## Parameters

- `k` : values of the wave numbers (one array for each cosines)
- `α` : strength of perturbations
- `σ` : variance of the Gaussian (one velocity vector for each gaussian).
- `μ` : mean value of the Gaussian (one velocity vector for each gaussian).
- `δ` : portion of each Gaussian 

# Example

```math
f(x,v_1,v_2)=\\frac{1}{2\\pi\\sigma_1\\sigma_2} \\exp \\Big( - \\frac{1}{2}
\\big( \\frac{v_1^2}{\\sigma_1^2} + \\frac{v_2^2}{\\sigma_2^2} \\big) 
\\Big) ( 1 + \\alpha \\cos(kx)),
```

```julia
df = CosSumGaussian{1,1,3}([[k]],[α], [[σ₁,σ₂]], [[μ₁,μ₂]])
```

"""
struct CosSumGaussian{D, V, S} <: AbstractCosGaussian

    dims        :: Tuple{Int64, Int64, Int64}
    params      :: CosGaussianParams

    function CosSumGaussian{D, V, S}( k :: Array{Vector{Float64}, 1}, 
                                   α :: Vector{Float64}, 
                                   σ :: Array{Vector{Float64}, 1}, 
                                   μ :: Array{Vector{Float64}, 1},
                                   δ :: Vector{Float64} = [1.0] ) where {D,V,S}

        dims   = (D, V, S)
        params = CosGaussianParams( dims, k, α, σ, μ, δ )

        new( dims, params )
    end

end

export SumCosGaussian

"""
    SumCosGaussian( dims, n_cos, n_gaussians, k, α, σ, μ, δ )
Data type for parameters of initial distribution
```math
(1+ \\sum_{i=1}^{n_{cos}} \\alpha_i \\cos(  k_i \\mathbf{x}))
\\cdot
\\sum_{j=1}^{n_{gaussians}} 
\\delta_j \\exp 
\\big( -\\frac{1}{2} 
\\frac{(\\mathbf{v}-\\mu_j)^2}{\\sigma_j^2} \\big)
```
## Parameters
- `k` : values of the wave numbers (Array of vectors for multiple cosines)
- `α` : strength of perturbations
- `σ` : variance of the Gaussian ( Array of vectors for multiple Gaussians)
- `μ` : mean value of the Gaussian ( Array multiple Gaussians)
- `normal` : Normalization constant of each Gaussian
- `n_gaussians` : Number of Gaussians
- `n_cos` : Number of cosines
- `δ` : portion of each Gaussian 

# Example

```math
f(x,v_1,v_2) = \\frac{1}{2\\pi\\sigma_1\\sigma_2} 
\\exp \\Big( - \\frac{1}{2} \\big( \\frac{v_1^2}{\\sigma_1^2}
 + \\frac{v_2^2}{\\sigma_2^2} \\big) \\Big) 
( 1 + \\alpha_1 \\cos(k_1 x) + \\alpha_2 \\cos(k_2 x) ),
```
```julia
df = SumCosGaussian{1,2}([[k₁],[k₂]], [α₁, α₂], [[σ₁,σ₂]], [[0.0,0.0]])

```
"""
struct SumCosGaussian{D,V,S} <: AbstractCosGaussian
    dims   :: Tuple{Int64,Int64,Int64}
    params :: CosGaussianParams

    function SumCosGaussian{D,V,S}( k :: Array{Array{Float64,1},1}, 
                                  α :: Vector{Float64}, 
                                  σ :: Array{Array{Float64,1},1}, 
                                  μ :: Array{Array{Float64,1},1},
                                  δ :: Array{Float64,1} = [1.0]
                                  ) where {D,V,S}

        dims   = (D, V, S)
        params = CosGaussianParams( dims, k, α, σ, μ, δ)
        new( dims, params )
    end

end

function eval_x_density( f :: CosSumGaussian, x :: Union{Float64,Vector{Float64}} )
    
    fval = 1.0
    for j=1:f.params.n_cos
       fval += f.params.α[j] * cos( sum(f.params.k[j] .* x) )
    end
    fval

end

"""
    eval_x_density( f, x )

evaluate the cosine part of the distribution function
"""
function eval_x_density( f :: SumCosGaussian, 
                         x :: Union{Float64,Vector{Float64}} )
    
    fval = 1.0
    for j=1:f.params.n_cos
       fval += f.params.α[j] * cos( sum(f.params.k[j] .* x) )
    end
    fval

end

"""
    eval_v_density( f, v )

evaluate the normal part of the distribution function
"""
function eval_v_density( f :: AbstractCosGaussian, 
                         v :: Union{Float64,Vector{Float64}} ) 

    fval = 0.0
    for j=1:f.params.n_gaussians
       fval += f.params.normal[j] * f.params.δ[j] .* exp( - 0.5 * 
               sum( ((v .- f.params.μ[j]) ./ f.params.σ[j]).^2))
    end
    fval

end 

function( f :: CosSumGaussian )( x, v )

    eval_x_density( f, x) * eval_v_density( f, v)

end

function( f :: SumCosGaussian )( x, v )

    eval_x_density( f, x) * eval_v_density( f, v)

end
