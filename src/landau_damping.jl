using Sobol

export LandauDamping, sample!

struct LandauDamping 

    alpha :: Float64
    kx    :: Float64

end

"""
    sample!( LandauDamping(α, kx), pg)

Sampling from a probability distribution to initialize
a Landau damping

```math
f_0(x,v,t) = \\frac{n_0}{2π v_{th}^2} ( 1 + \\alpha cos(k_x x))
 exp( - \\frac{v_x^2+v_y^2}{2 v_{th}^2})
```
The newton function solves the equation ``P(x)-r=0`` with Newton’s method
```math
    x^{n+1} = x^n – (P(x)-(2\\pi r / k)/f(x) 
```
with 
```math
P(x) = \\int_0^x (1 + \\alpha cos(k_x y)) dy = x + \\frac{\\alpha}{k_x} sin(k_x x)
```

"""
function sample!( d :: LandauDamping, pg :: ParticleGroup{1,2} )

    alpha, kx = d.alpha, d.kx
    
    function newton(r)
        x0, x1 = 0.0, 1.0
        r *= 2π / kx
        while (abs(x1-x0) > 1e-12)
            p = x0 + alpha * sin( kx * x0) / kx 
            f = 1 + alpha * cos( kx * x0)
            x0, x1 = x1, x0 - (p - r) / f
        end
        x1
    end
    
    s = SobolSeq(2)
    nbpart = pg.n_particles

    for i=1:nbpart
        v = sqrt(-2 * log( (i-0.5)/nbpart))
        r1, r2 = Sobol.next!(s)
        θ = r1 * 2π
        set_x(pg,  i, newton(r2))
        set_v(pg,  i, [v * cos(θ), v * sin(θ)])
        set_weights( pg, i, 2*pi/kx/nbpart)
    end


end
