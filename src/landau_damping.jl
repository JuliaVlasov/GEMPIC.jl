import Sobol


"""
```math
    f(x) = \\frac{1]{2\\pi} (1 + \\alpha cos(k x)) 
    \\exp \\big( - \\frac{ (v - \\mu)^2 }{ \\sigma^2} \\big)
```
"""
struct LandauDamping

    α :: Float64
    k :: Float64
    σ :: Float64

end

"""
Input r is a random number ``\\in [0,1]``

```math
    f(x) = 1 + \\alpha cos(k x)
```
on some domain ``[0, 2\\pi/k]``

Solve the equation ``P(x)-r=0`` with Newton’s method

```math
    x^{n+1} = x^n – (P(x)-(2\\pi r / k)/f(x) 
```

with 
```math
P(x) = \\int_0^x (1 + \\alpha cos(k y)) dy
```
```math
P(x) = x + \\frac{\\alpha}{k} sin (k x)
```
"""
function newton( ld :: LandauDamping, r)
    x0, x1 = 0.0, 1.0
    alpha, k = ld.α, ld.k
    r *= 2π / k
    while (abs(x1-x0) > 1e-12)
        p = x0 + alpha * sin( k * x0) / k 
        f = 1 + alpha * cos( k * x0)
        x0, x1 = x1, x0 - (p - r) / f
    end
    x1
end

"""
    sample!( landau, xp, yp, domain )

Particle sampling for Landau damping initial distribution
function. 
- `xp` :: preallocated array containing particles positions
- `vp` :: preallocated array containing particles velocities
"""
function sample!( ld :: LandauDamping, 
                  xp :: Vector{Float64}, 
                  vp :: Vector{Float64},
                  domain :: Vector{Float64} )
    
   nbpart = length(xp)
   @assert length(xp) == length(vp)
    
   s = SobolSeq(2)

   for k=0:nbpart-1

	  v = ld.σ * sqrt(-2 * log( (ld.k+0.5)/nbpart))
      r1, r2 = Sobol.next!(s)
      θ = r1 * 2π
      xp[k+1] =  newton(r2)
      vp[k+1] =  v * sin(θ)

   end

end

#xp, vp = landau(100000);
# -
#
#p = histogram([xp,vp], normalize=true, bins = 100,  layout=(2,1), lab = "draws")
#plot!(p[1,1], x-> (1+0.1*cos(0.5*x))/4π, 0., 4π)
#plot!(p[2,1], x-> (exp(-x^2/2))/sqrt(2π), -6, 6)
