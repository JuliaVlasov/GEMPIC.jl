abstract type CosGaussianDistributions end

export CosSumOneGaussian

"""
    CosGaussian( dims, n_cos, n_gaussians, kx, alpha, v_thermal, v_mean, δ )

Data type for parameters of initial distribution

```math
(1+ \\cos( \\sum_{n_{cos}} kx_i * x_i)) * \\exp \\big( -\\frac{1}{2} \\sum_{n_{gaussians} \\frac{(v-v_{mean})^2}{v_{thermal}^2} \\big)
```

## Parameters

- `kx`          : values of the wave numbers (first index dimension, second index for multiple cosines)
- `alpha`       : strength of perturbations
- `v_thermal`   : variance of the Gaussian ( first index velocity dimension, second index multiple Gaussians)
- `v_mean`      : mean value of the Gaussian ( first index velocity dimension, second index multiple Gaussians)
- `delta`       : Portion of each Gaussian
- `normal`      : Normalization constant of each Gaussian
- `n_gaussians` : Number of Gaussians
- `n_cos`       : Number of cosines

# Example

```math
f(x,i\\mathbf{v},t=0)=\\frac{1}{2\\pi\\sigma_1\\sigma_2} \\exp \\Big( - \\frac{1}{2} \\big( \\frac{v_1^2}{\\sigma_1^2} + \\frac{v_2^2}{\\sigma_2^2} \\big) \\Big) ( 1 + \\alpha \\cos(kx)), \\qquad x \\in [0, 2\\pi / k),
```

```julia
kx        = hcat([1.25])
v_thermal = hcat([σ₁, σ₂])
v_mean    = hcat([0.0, 0.0])
delta     = [1.0]

```


"""
struct CosSumOneGaussian

    dims        :: Tuple{Int64, Int64}
    n_cos       :: Int64
    n_gaussians :: Int64
    kx          :: Array{Float64, 2} 
    alpha       :: Vector{Float64}
    v_thermal   :: Array{Float64, 2}
    v_mean      :: Array{Float64, 2}
    delta       :: Vector{Float64}
    normal      :: Vector{Float64}

    function CosSumOneGaussian( dims        :: Tuple{Int64, Int64},
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

        new( dims, n_cos, n_gaussians, kx, alpha, v_thermal, 
             v_mean, delta, normal )
    end

end

export eval_x_density

function eval_x_density( self :: CosSumOneGaussian, x )
    
    fval = 1.0
    for j=1:self.n_cos
       fval += self.alpha[j] * cos( sum(self.kx[:,j] .* x) )
    end
    fval

end
  
export eval_v_density

function eval_v_density( self :: CosSumOneGaussian, v ) 

    fval = 0.0
    for j=1:self.n_gaussians
       fval += self.normal[j] * self.delta[j] * exp( - 0.5 * 
               sum( ((v .- self.v_mean[:,j]) ./ self.v_thermal[:,j]).^2))
    end
    fval

end 

export eval_xv_density
function eval_xv_density( self :: CosSumOneGaussian, x, v ) 

    eval_x_density(self, x) * eval_v_density(self, v)

end 

function (self :: CosSumOneGaussian)( x :: Float64, v :: Float64 ) 

    eval_x_density( self, x) * eval_v_density( self, v)

end

function (self :: CosSumOneGaussian)( x::Vector{Float64}, v::Vector{Float64} ) 

    f = zeros(Float64, (length(x), length(v)))

    for j in eachindex(v), i in eachindex(x)
       f[i,j] = eval_x_density( self, x[i]) * eval_v_density( self, v[j])
    end

    f
end
