@testset " spline pp 1d " begin

import GEMPIC: SplinePP, b_to_pp, uniform_bsplines_eval_basis, horner_1d

ncells = 8
degree = 3
  
domain    = [0., 2π]
delta_x   = (domain[2] - domain[1]) / ncells
 
b_coeffs  = [8.9559045978573115E-002,  0.53453935506698202, 
             0.12723516577130212    ,  0.12700877094628704,
             0.19068356307822054    ,  0.23912829562467719,
             0.64734047946930895    ,  0.46282228199849473]

xp = 4.8141437173987800
 
spline_pp = SplinePP( degree, ncells)

pp_coeffs = b_to_pp(spline_pp, ncells, b_coeffs)

xi    = (xp - domain[1])/delta_x
index = floor(Int64, xi)+1
xi    = xi - (index-1)
res   = horner_1d(degree, pp_coeffs, xi, index)
    
index = index - degree
    
val = uniform_bsplines_eval_basis(degree, xi)

res2 = 0.0

for i = 1:degree+1 
   index1d = (index+i-2) % ncells + 1
   res2    += b_coeffs[index1d] * val[i]
end 

@test abs(res-res2) < 1e-15

@test spline_pp.degree - degree < 1e-15

xp    = rand()
index = 1
res = horner_1d(degree, pp_coeffs, xp, index)

res2 = 0.
for i=1:degree+1
    res2 += pp_coeffs[i,1] * xp^((degree+1)-i)
end

@test abs(res-res2) < 1e-12
     
end
