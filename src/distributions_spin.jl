export SpinCosSumGaussian

"""
    SpinCosGaussianParams( dims, k, α, σ, μ, δ ) 

Parameters of a distribution with is a product of a Cosine 
distribution along x and a Normal distribution along v.

- `n_gaussians` : Number of Gaussians
- `n_cos`       : Number of cosines
- `normal`      : Normalization constant of each Gaussian
"""
struct SpinCosGaussianParams

    dims        :: Tuple{Int64, Int64, Int64}
    n_cos       :: Int64
    n_gaussians :: Int64
    k           :: Array{Vector{Float64}, 1}
    α           :: Vector{Float64}
    σ           :: Array{Vector{Float64}, 1}
    μ           :: Array{Vector{Float64}, 1}
    normal      :: Vector{Float64}
    δ           :: Vector{Float64}

    function SpinCosGaussianParams( dims :: Tuple{Int64, Int64, Int64},
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

export SpinCosSumGaussian

"""
    SpinCosSumGaussian{D,V}( n_cos, n_gaussians, k, α, σ, μ, δ )

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
struct SpinCosSumGaussian{D, V, S} <: AbstractCosGaussian

    dims        :: Tuple{Int64, Int64, Int64}
    params      :: SpinCosGaussianParams

    function SpinCosSumGaussian{D, V, S}( k :: Array{Vector{Float64}, 1}, 
                                   α :: Vector{Float64}, 
                                   σ :: Array{Vector{Float64}, 1}, 
                                   μ :: Array{Vector{Float64}, 1},
                                   δ :: Vector{Float64} = [1.0] ) where {D,V,S}

        dims   = (D, V, S)
        params = SpinCosGaussianParams( dims, k, α, σ, μ, δ )

        new( dims, params )
    end

end

function eval_x_density( f :: SpinCosSumGaussian, x :: Union{Float64,Vector{Float64}} )
    
    fval = 1.0
    for j=1:f.params.n_cos
       fval += f.params.α[j] * cos( sum(f.params.k[j] .* x) )
    end
    fval

end

function( f :: SpinCosSumGaussian )( x, v )

    eval_x_density( f, x) * eval_v_density( f, v)

end
