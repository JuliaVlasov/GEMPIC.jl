using StaticArrays

const  inv_2   = 1. / 2.
const  inv_3   = 1. / 3.
const  inv_4   = 1. / 4.
const  inv_6   = 1. / 6.
const  inv_8   = 1. / 8.
const  inv_12  = 1. / 12.
const  inv_18  = 1. / 18.
const  inv_20  = 1. / 20.
const  inv_24  = 1. / 24.
const  inv_30  = 1. / 30.
const  inv_36  = 1. / 36.
const  inv_48  = 1. / 48.
const  inv_72  = 1. / 72.
const  inv_120 = 1. / 120.
const  inv_144 = 1. / 144.
const  inv_720 = 1. / 720.

"""
    SplinePP( degree, ncells)

- `degree` : degree of 1d spline
- `poly_coeffs` : `poly_coeffs[i,j]` coefficient of ``x^{deg+1-j}`` for ith B-spline function  size= (degree+1, degree+1)
- `poly_coeffs_fp` : `poly_coeffs[i,j]` coefficient of ``x^{deg+1-j}`` for ith B-spline function  size= (degree+1, degree+1)
- `ncells` : number of gridcells
- `scratch_b` : scratch data for b_to_pp-converting
- `scratch_p` : scratch data for b_to_pp-converting
"""
mutable struct SplinePP

    degree         :: Int64
    poly_coeffs    :: Array{Float64,2}
    poly_coeffs_fp :: Array{Float64,2}
    ncells         :: Int64
    scratch_b      :: Vector{Float64}
    scratch_p      :: Vector{Float64}

    function SplinePP( degree, ncells)
    
         @assert (ncells >= degree )
         poly_coeffs    = zeros(Float64, (degree+1,degree+1))
         poly_coeffs_fp = zeros(Float64, (degree+1,degree+1))
         scratch_b      = zeros(Float64, (degree+1))
         scratch_p      = zeros(Float64, (degree+1))
    
    if degree == 1

         poly_coeffs    = SMatrix{2,2}( -1., 1., 1., 0. ) 
         poly_coeffs_fp = SMatrix{2,2}( -inv_2, 1., inv_2, 0. )
       
    elseif degree == 2

         poly_coeffs    = SMatrix{3,3}(inv_2, -1., inv_2 , -1., 1., 
                                   inv_2, inv_2, 0., 0.)
         poly_coeffs_fp = SMatrix{3,3}(inv_6, -inv_2, inv_2 , -inv_3, 
                                   inv_2, inv_2, inv_6, 0., 0)

    elseif degree == 3

        poly_coeffs = SMatrix{4,4}( -inv_6, inv_2, -inv_2, inv_6
                              ,  inv_2,   -1.,     0., 4*inv_6
                              , -inv_2, inv_2,  inv_2, inv_6
                              ,  inv_6,    0.,     0., 0.)
       
        poly_coeffs_fp = SMatrix{4,4}(- inv_24, inv_6, -inv_4, inv_6
                                 ,  inv_8, -inv_3,     0., 4*inv_6
                                 , -inv_8,  inv_6,  inv_4, inv_6
                                 ,  inv_24,    0.,     0., 0.)

    elseif degree == 4

        poly_coeffs = SMatrix{5,5}(inv_24,-inv_6, inv_4,-inv_6, inv_24
                             ,- inv_6, inv_2,-inv_4,-inv_2, 11*inv_24
                             ,  inv_4,-inv_2,-inv_4, inv_2, 11*inv_24
                             ,- inv_6, inv_6, inv_4, inv_6, inv_24
                             , inv_24,    0.,    0.,    0., 0.   )

        poly_coeffs_fp = SMatrix{5,5}( inv_120,- inv_24, inv_12,-inv_12, inv_24
                                , - inv_30,  inv_8,-inv_12,-inv_2, 11*inv_24
                                ,   inv_20,- inv_8,-inv_12,inv_4,11*inv_24
                                , - inv_30,  inv_24,inv_12,inv_12,inv_24
                                ,   inv_120,     0.,    0.,    0., 0.)

    elseif degree == 5

        poly_coeffs = SMatrix{6,6}(-inv_120,inv_24,-inv_12,inv_12,-inv_24,inv_120
                               ,inv_24,-inv_6,inv_6,inv_6,-5*inv_12, 26*inv_120 
                               ,-inv_12,inv_4,0.,-inv_2,0.,11*inv_20
                               ,inv_12,-inv_6,-inv_6,inv_6,5*inv_12,26*inv_120 
                               ,-inv_24,inv_24,inv_12,inv_12,inv_24,inv_120
                               ,inv_120,0.,0.,0.,0.,0.)

        poly_coeffs_fp = SMatrix{6,6}(-inv_720,inv_120,-inv_48,inv_36,-inv_48,inv_120
                                  , inv_144,-inv_30,inv_24,inv_18,-5*inv_24, 26*inv_120
                                  ,-inv_72,inv_20,0.,-inv_6,0.,11*inv_20
                                  ,inv_72,-inv_30,-inv_24,inv_18,5*inv_24,26*inv_120
                                  ,-inv_144,inv_120,inv_48,inv_36,inv_48,inv_120
                                  ,inv_720,0.,0.,0.,0.,0.) 
    else

       throw(ArgumentError(" degree $degree not implemented"))

    end

    new( degree, poly_coeffs, poly_coeffs_fp, ncells, scratch_b, scratch_p)

  end 
     
end 

"""
    b_to_pp( SplinePP, ncells, b_coeffs)

Convert 1d spline in B form to spline in pp form with 
periodic boundary conditions
"""
function b_to_pp( self :: SplinePP, ncells :: Int64, b_coeffs :: Vector{Float64})

    degp1     = self.degree+1
    pp_coeffs = zeros(Float64, (degp1,ncells)) 
    coeffs    = zeros(Float64, degp1) 
       
    for i=1:self.degree   
        coeffs .= vcat(b_coeffs[end-self.degree+i:end],b_coeffs[1:i]) 
        for j=1:degp1
            pp_coeffs[j, i] = sum(coeffs .* self.poly_coeffs[j,:])
        end
    end
    
    for i=self.degree+1:ncells
        coeffs .= b_coeffs[i-self.degree:i]
        for j=1:degp1
            pp_coeffs[j, i] = sum(coeffs .* self.poly_coeffs[j,:])
        end
    end

    pp_coeffs

end

"""
    horner_1d(degree, pp_coeffs, x, index)

Perform a 1d hornerschema on the pp_coeffs at index
"""
function horner_1d(degree :: Int, pp_coeffs, x :: Float64, index :: Int)
    
    res = pp_coeffs[1,index]
    for i=1:degree
       res = res * x + pp_coeffs[i+1,index]
    end
    res

end

"""
    horner_primitive_1d(val, degree, pp_coeffs, x)

Perform a 1d hornerschema on the pp_coeffs evaluate at x
"""
function horner_primitive_1d(val :: Vector{Float64}, degree, pp_coeffs, x)

  for i in eachindex(val)

     val[i] = horner_1d(degree, pp_coeffs, x, i) * x

  end

end 

