"""
    uniform_bsplines_eval_basis( spline_degree, normalized_offset, bspl )
  
# UNIFORM B-SPLINE FUNCTIONS

## Evaluate all non vanishing uniform B-Splines in unit cell. 

Returns an array with the values of the b-splines of the 
requested degree, evaluated at a given cell offset. The cell size is
normalized between 0 and 1, thus the offset given must be a number
between 0 and 1.

Output: 

```math
bspl(1:d+1)= B_d(-(d+1)/2+d+x),...,B_d(-(d+1)/2+x)
```
 
with ``d``=`spline_degree` and ``x``=`normalized_offset`
where ``B_d=B_{d-1}*B_0`` and ``B_0=1_[-1/2,1/2]`` and `*` is convolution
the following FORTRAN code can be used for comparison with 
[deboor](http://pages.cs.wisc.edu/~deboor/)

```fortran
do i=-d,d+1
    t(i+d+1)=real(i,8)
end do
call bsplvb(t,d+1,1,normalized_offset,d+1,out)
```

We also have the property (from the symmetry of the B-spline)
```math
out[1:d+1]= B_d(-(d+1)/2+xx),...,B_d(-(d+1)/2+d+xx),..., 
```
where ``xx=1-`` `normalized_offset`

"""
function uniform_bsplines_eval_basis(spline_degree::Int, normalized_offset::Float64)

    @assert spline_degree >= 0
    @assert normalized_offset >= 0.0
    @assert normalized_offset <= 1.0

    bspl = zeros(Float64, spline_degree + 1)

    bspl[1] = 1.0
    @inbounds for j = 1:spline_degree
        xx = -normalized_offset
        j_real = Float64(j)
        inv_j = 1.0 / j_real
        saved = 0.0
        for r = 0:j-1
            xx = xx + 1.0
            temp = bspl[r+1] * inv_j
            bspl[r+1] = saved + xx * temp
            saved = (j_real - xx) * temp
        end
        bspl[j+1] = saved
    end

    bspl

end

function uniform_bsplines_eval_basis!(
    bspl::Vector{Float64},
    spline_degree::Int,
    normalized_offset::Float64,
)
    bspl[1] = 1.0
    @inbounds for j = 1:spline_degree
        xx = -normalized_offset
        j_real = Float64(j)
        inv_j = 1.0 / j_real
        saved = 0.0
        for r = 0:j-1
            xx = xx + 1.0
            temp = bspl[r+1] * inv_j
            bspl[r+1] = saved + xx * temp
            saved = (j_real - xx) * temp
        end
        bspl[j+1] = saved
    end

end

function uniform_bsplines_eval_basis!(
    bspl::Array{Float64,2},
    spline_degree::Int,
    normalized_offset_x::Float64,
    normalized_offset_y::Float64,
)
    bspl[1, 1] = 1.0
    bspl[1, 2] = 1.0
    @inbounds for j = 1:spline_degree
        xx = -normalized_offset_x
        yy = -normalized_offset_y
        j_real = Float64(j)
        inv_j = 1.0 / j_real
        saved_x = 0.0
        saved_y = 0.0
        for r = 0:j-1
            xx = xx + 1
            yy = yy + 1
            temp_x = bspl[r+1, 1] * inv_j
            temp_y = bspl[r+1, 2] * inv_j
            bspl[r+1, 1] = saved_x + xx * temp_x
            bspl[r+1, 2] = saved_y + yy * temp_y
            saved_x = (j_real - xx) * temp_x
            saved_y = (j_real - yy) * temp_y
        end
        bspl[j+1, 1] = saved_x
        bspl[j+1, 2] = saved_y
    end

end

export eval_uniform_periodic_spline_curve

"""
    eval_uniform_periodic_spline_curve( degree, scoef )

Evaluate uniform periodic spline curve defined by coefficients scoef at 
knots (which are the grid points) 

"""
function eval_uniform_periodic_spline_curve(degree::Int, scoef::Vector{Float64})

    # get bspline values at knots
    bspl = uniform_bsplines_eval_basis(degree, 0.0)
    n = length(scoef)
    sval = similar(scoef)
    for i = 1:n
        val = 0.0
        for j = 1:degree
            imj = mod1(i - j + n, n)
            val = val + bspl[j] * scoef[imj]
        end
        sval[i] = val
    end

    sval

end
