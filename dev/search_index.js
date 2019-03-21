var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Documentation",
    "title": "Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "#GEMPIC.Maxwell1DFEM",
    "page": "Documentation",
    "title": "GEMPIC.Maxwell1DFEM",
    "category": "type",
    "text": "maxwell_solver = MaxwellFEM1D( domain, ncells, degree )\n\n1D Maxwell spline finite element solver on a periodic grid\n\nLx                   : length of Periodic domain\ndelta_x              : cell size\nn_dofs               : number of cells (and grid points)\nsdeg0              : spline degree 0-forms\nsdeg1              : spline degree 1-forms\nmass_0               : coefficients of 0-form mass matrix\nmass_1               : coefficients of 1-form mass matrix\neig_mass0            : eigenvalues of circulant 0-form mass matrix\neig_mass1            : eigenvalues of circulant 1-form mass matrix\neigweakampere      : eigenvalues of circulant update matrix for Ampere\neigweakpoisson     : eigenvalues of circulant update matrix for Poisson\nplan_fw              : fft plan (forward)\nplan_bw              : fft plan (backward)\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.ParticleMeshCoupling",
    "page": "Documentation",
    "title": "GEMPIC.ParticleMeshCoupling",
    "category": "type",
    "text": "Kernel smoother for 2d with splines of arbitrary degree placed on a uniform mesh. Spline with index i starts at point i\n\nValue of grid spacing along both directions.\nDefinition of the domain: domain(1,1:2) = x1min, x1max\nNumber of particles of underlying PIC method (processor local)\nDegree of smoothing kernel spline\nNumber of intervals where spline non zero (spline_degree + 1)\nScaling factor depending on whether Galerkin or collocation\nNumber of quadrature points\nscratch data for spline evaluation\nmore scratch data for spline evaluation\nquadrature weights and points\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.compute_b_from_e-NTuple{4,Any}",
    "page": "Documentation",
    "title": "GEMPIC.compute_b_from_e",
    "category": "method",
    "text": "Compute Bz from Ey using strong 1D Faraday equation for spline coefficients\n\nB_z^new(x_j) = B_z^old(x_j) - fracDelta tDelta x (E_y(x_j) - E_y(x_j-1)\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.compute_e_from_b-NTuple{4,Any}",
    "page": "Documentation",
    "title": "GEMPIC.compute_e_from_b",
    "category": "method",
    "text": "compute Ey from Bz using weak Ampere formulation \n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.compute_e_from_j-NTuple{4,Any}",
    "page": "Documentation",
    "title": "GEMPIC.compute_e_from_j",
    "category": "method",
    "text": "Compute E_i from j_i integrated over the time interval using weak Ampere formulation\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.compute_rhs_from_function-NTuple{4,Any}",
    "page": "Documentation",
    "title": "GEMPIC.compute_rhs_from_function",
    "category": "method",
    "text": "computerhsfromfunction(self, func, degree, coefsdofs)\n\nCompute the FEM right-hand-side for a given function f and periodic splines of given degree.\n\nIts components are int f N_i dx where N_i is the B-spline starting at x_i. \n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.eval_uniform_periodic_spline_curve-Tuple{Int64,Array{Float64,1}}",
    "page": "Documentation",
    "title": "GEMPIC.eval_uniform_periodic_spline_curve",
    "category": "method",
    "text": "evaluniformperiodicsplinecurve( degree, scoef )\n\nEvaluate uniform periodic spline curve defined by coefficients scoef at  knots (which are the grid points) \n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.l2norm_squared-Tuple{Any,Any,Any}",
    "page": "Documentation",
    "title": "GEMPIC.l2norm_squared",
    "category": "method",
    "text": "Compute square of the L2norm\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.l2projection-NTuple{4,Any}",
    "page": "Documentation",
    "title": "GEMPIC.l2projection",
    "category": "method",
    "text": "Compute the L2 projection of a given function f on periodic splines of given degree\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.SplinePP",
    "page": "Documentation",
    "title": "GEMPIC.SplinePP",
    "category": "type",
    "text": "degree         : degree of 1d spline\npolycoeffs    : polycoeffs[i,j] coefficient of x^{deg+1-j} for ith B-spline function  size= (degree+1, degree+1)\npolycoeffsfp : poly_coeffs[i,j] coefficient of x^{deg+1-j} for ith B-spline function  size= (degree+1, degree+1)\nncells        : number of gridcells\nscratchb      : scratch data for bto_pp-converting\nscratchp      : scratch data for bto_pp-converting\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.b_to_pp-Tuple{Any,Any,Any}",
    "page": "Documentation",
    "title": "GEMPIC.b_to_pp",
    "category": "method",
    "text": "Convert 1d spline in B form to spline in pp form with periodic boundary conditions\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.get_charge-Tuple{ParticleGroup1D2V,Any}",
    "page": "Documentation",
    "title": "GEMPIC.get_charge",
    "category": "method",
    "text": "Get charge of particle (q * particle_weight)\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.get_mass-Tuple{ParticleGroup1D2V,Any}",
    "page": "Documentation",
    "title": "GEMPIC.get_mass",
    "category": "method",
    "text": "Get mass of particle (m * particle_weight)\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.get_v-Tuple{ParticleGroup1D2V,Any}",
    "page": "Documentation",
    "title": "GEMPIC.get_v",
    "category": "method",
    "text": "Get velocities\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.get_weights-Tuple{ParticleGroup1D2V,Any}",
    "page": "Documentation",
    "title": "GEMPIC.get_weights",
    "category": "method",
    "text": "Get particle weights\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.get_x-Tuple{ParticleGroup1D2V,Any}",
    "page": "Documentation",
    "title": "GEMPIC.get_x",
    "category": "method",
    "text": "Get position\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.horner_1d-Tuple{Int64,Any,Float64,Int64}",
    "page": "Documentation",
    "title": "GEMPIC.horner_1d",
    "category": "method",
    "text": "Perform a 1d hornerschema on the pp_coeffs at index\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.set_common_weight-Tuple{ParticleGroup1D2V,Any}",
    "page": "Documentation",
    "title": "GEMPIC.set_common_weight",
    "category": "method",
    "text": "Set the common weight\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.set_v-Tuple{ParticleGroup1D2V,Any,Any}",
    "page": "Documentation",
    "title": "GEMPIC.set_v",
    "category": "method",
    "text": "Set velocity of particle @ i\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.set_weights-Tuple{ParticleGroup1D2V,Any,Any}",
    "page": "Documentation",
    "title": "GEMPIC.set_weights",
    "category": "method",
    "text": "Set weights of particle @ i\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.set_x-Tuple{ParticleGroup1D2V,Any,Any}",
    "page": "Documentation",
    "title": "GEMPIC.set_x",
    "category": "method",
    "text": "Set position of particle @ i\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.uniform_bsplines_eval_basis-Tuple{Int64,Float64}",
    "page": "Documentation",
    "title": "GEMPIC.uniform_bsplines_eval_basis",
    "category": "method",
    "text": "uniform_bsplines_eval_basis( spline_degree, normalized_offset, bspl )\n\nUNIFORM B-SPLINE FUNCTIONS\n\nEvaluate all non vanishing uniform B-Splines in unit cell.\n\nReturns an array with the values of the b-splines of the  requested degree, evaluated at a given cell offset. The cell size is normalized between 0 and 1, thus the offset given must be a number between 0 and 1.\n\nOutput: bspl(1:d+1)= Bd(-(d+1)/2+d+x),...,Bd(-(d+1)/2+x)  with d=splinedegree and x=normalizedoffset where Bd=B{d-1}*B0 and B0=1_[-1/2,1/2] and * is convolution the following code can be used for comparison with deboor\n\ndo i=-d,d+1\n    t(i+d+1)=real(i,8)\nend do\ncall bsplvb(t,d+1,1,normalized_offset,d+1,out)\n\nWe also have the property (from the symmetry of the B-spline) out(1:d+1)= Bd(-(d+1)/2+xx),...,Bd(-(d+1)/2+d+xx),...,  where xx=1-normalized_offset\n\n\n\n\n\n"
},

{
    "location": "#GEMPIC.jl-Documentation-1",
    "page": "Documentation",
    "title": "GEMPIC.jl Documentation",
    "category": "section",
    "text": "Modules = [GEMPIC]\nOrder   = [:type, :function]"
},

{
    "location": "contents/#",
    "page": "Contents",
    "title": "Contents",
    "category": "page",
    "text": ""
},

{
    "location": "contents/#Contents-1",
    "page": "Contents",
    "title": "Contents",
    "category": "section",
    "text": ""
},

{
    "location": "contents/#Index-1",
    "page": "Contents",
    "title": "Index",
    "category": "section",
    "text": ""
},

]}
