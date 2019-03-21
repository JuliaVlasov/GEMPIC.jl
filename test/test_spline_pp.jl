@testset " spline pp 1d " begin

import GEMPIC: SplinePP
using Random

ncells = 8
degree = 3
  
b_coeffs  = zeros(Float64, ncells)
pp_coeffs = zeros(Float64,(degree+1,ncells))
val       = zeros(Float64, degree+1)

domain    = [0., 2π]
delta_x   = (domain[2] - domain[1]) / ncells
 
Random.seed!(42)
b_coeffs  = rand(Float64,ncells)
 
spline_pp = SplinePP( degree, ncells)

b_to_pp(spline_pp, n_cells, b_coeffs, pp_coeffs)

xp = rand() * (domain[2]-domain[1])

xi    = (xp - domain[1])/delta_x
index = floor(xi)+1
xi    = xi - (index-1)
res   = horner_1d(degree, pp_coeffs, xi, index)
    
index = index - degree
    
val = uniform_bsplines_eval_basis(degree, xi)
  
res2 = 0.0

for i = 1:degree+1 
   index1d = modulo(index+i-2, n_cells)+1
   res2 = res2 + b_coeffs(index1d) * val(i)
end 

@test abs(res-res2) < 1e-15

@test spline_pp.degree - degree < 1e-15

res = horner_1d(degree, pp_coeffs, rand(), 1)

res2 = 0.
for i=1:degree+1
   res2 += pp_coeffs[i,1] * xp^((degree+1)-i)
end

@test abs(res-res2) < 1e-12
     
end
