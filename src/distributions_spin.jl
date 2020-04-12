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
struct SpinCosSumGaussian <: AbstractCosGaussian

    dims        :: Tuple{Int64, Int64, Int64}
    params      :: CosGaussianParams

    function SpinCosSumGaussian( k :: Array{Vector{Float64}, 1}, 
                                 α :: Vector{Float64}, 
                                 σ :: Array{Vector{Float64}, 1}, 
                                 μ :: Array{Vector{Float64}, 1},
                                 δ :: Vector{Float64} = [1.0] )

        dims   = (1, 1)
        params = CosGaussianParams( dims, k, α, σ, μ, δ )

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
