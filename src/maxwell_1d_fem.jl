using LinearAlgebra
using FFTW
using FastGaussQuadrature

abstract type AbstractMaxwellSolver end

export Maxwell1DFEM

"""
    maxwell_solver = MaxwellFEM1D( domain, ncells, degree )

1D Maxwell spline finite element solver on a periodic grid

- `Lx`                   : length of Periodic domain
- `delta_x`              : cell size
- `n_dofs`               : number of cells (and grid points)
- `s_deg_0`              : spline degree 0-forms
- `s_deg_1`              : spline degree 1-forms
- `mass_0`               : coefficients of 0-form mass matrix
- `mass_1`               : coefficients of 1-form mass matrix
- `eig_mass0`            : eigenvalues of circulant 0-form mass matrix
- `eig_mass1`            : eigenvalues of circulant 1-form mass matrix
- `eig_weak_ampere`      : eigenvalues of circulant update matrix for Ampere
- `eig_weak_poisson`     : eigenvalues of circulant update matrix for Poisson
- `plan_fw`              : fft plan (forward)
- `plan_bw`              : fft plan (backward)

"""
mutable struct Maxwell1DFEM <: AbstractMaxwellSolver

    Lx               :: Float64  
    delta_x          :: Float64     
    n_dofs           :: Int32   
    s_deg_0          :: Int32
    s_deg_1          :: Int32
    mass_0           :: Vector{Float64}      
    mass_1           :: Vector{Float64}      
    eig_mass0        :: Vector{Float64}   
    eig_mass1        :: Vector{Float64}   
    eig_weak_ampere  :: Vector{Float64}  
    eig_weak_poisson :: Vector{Float64} 

    plan_fw :: FFTW.FFTWPlan
    plan_bw :: FFTW.FFTWPlan
    work    :: Vector{Float64}
    wsave   :: Vector{Float64}

    function Maxwell1DFEM( domain, ncells :: Int, degree :: Int )

        n_dofs  = ncells
        Lx      = domain[2] - domain[1]
        delta_x = Lx / n_dofs
        s_deg_0 = degree
        s_deg_1 = degree - 1

        mass_0  = zeros(Float64, s_deg_0+1)
        mass_1  = zeros(Float64, s_deg_0)

        if s_deg_0 == 1 # linear and constant splines

            # Upper diagonal coeficients  of linear spline mass matrix (from Eulerian numbers)
            mass_0[1] = 4.0/6.0
            mass_0[2] = 1.0/6.0
            # Upper diagonal coeficients  of constant spline mass matrix
            mass_1[1] = 1.0 

        elseif s_deg_0 == 2 # quadratic and linear splines

            # Upper diagonal coeficients  of quadratic spline mass matrix (from Eulerian numbers)
            mass_0[1] = 66.0/120.0
            mass_0[2] = 26.0/120.0
            mass_0[3] = 1.0/120.0
            # Upper diagonal coeficients  of linear spline mass matrix (from Eulerian numbers)
            mass_1[1] = 4.0/6.0 
            mass_1[2] = 1.0/6.0

        elseif s_deg_0 == 3

            # Upper diagonal coeficients  of cubic spline mass matrix (from Eulerian numbers)
            mass_0[1] = 2416.0/5040.0 
            mass_0[2] = 1191.0/5040.0
            mass_0[3] = 120.0/5040.0
            mass_0[4] = 1.0/5040.0
            # Upper diagonal coeficients  of quadratic spline mass matrix (from Eulerian numbers)
            mass_1[1] = 66.0/120.0 
            mass_1[2] = 26.0/120.0
            mass_1[3] = 1.0/120.0

        else

            throw( ArgumentError("Wrong value of degree = $degree  (1,2 or 3)") )

        end 

        eig_mass0         = zeros(Float64, n_dofs)
        eig_mass1         = zeros(Float64, n_dofs)
        eig_weak_ampere   = zeros(Float64, n_dofs)
        eig_weak_poisson  = zeros(Float64, n_dofs)

        work     = zeros(Float64, n_dofs)
        wsave    = zeros(Float64, n_dofs)
        plan_fw  = FFTW.plan_r2r(work,  FFTW.R2HC)
        plan_bw  = FFTW.plan_r2r(wsave, FFTW.HC2R)

        # Compute eigenvalues of circulant Ampere update matrix M_0^{-1} D^T M_1
        # and circulant Poisson Matrix (D^T M_1 D)^{-1}
        # zero mode vanishes due to derivative matrix D^T

        eig_weak_ampere[1]  = 0.0 
        eig_weak_poisson[1] = 0.0  # Matrix is not invertible: 0-mode is set to 0
        eig_mass0[1]        = 1.0  # sum of coefficents is one
        eig_mass1[1]        = 1.0  # sum of coefficents is one

        for k = 1:n_dofs÷2-1

            coef0 = mass_0[1]
            coef1 = mass_1[1]
            for j = 1:s_deg_0-1
               cos_mode = cos(2*pi*j*k/n_dofs)
               coef0 = coef0 + 2 * mass_0[j+1] * cos_mode
               coef1 = coef1 + 2 * mass_1[j+1] * cos_mode
            end
            # add last term for larger matrix
            j = s_deg_0
            coef0 = coef0 + 2 * mass_0[j+1]*cos(2*pi*j*k/n_dofs)
            # compute eigenvalues
            eig_mass0[k+1] = coef0 # real part
            eig_mass0[n_dofs-k+1] = 0.0 # imaginary part
            eig_mass1[k+1] = coef1 # real part
            eig_mass1[n_dofs-k+1] = 0.0 # imaginary part
            cos_mode = cos(2*pi*k/n_dofs)
            sin_mode = sin(2*pi*k/n_dofs)
            eig_weak_ampere[k+1] =  (coef1 / coef0) * (1-cos_mode) # real part
            eig_weak_ampere[n_dofs-k+1] =  -(coef1 / coef0) * sin_mode   # imaginary part
            eig_weak_poisson[k+1] = 1.0 / (coef1 * ((1-cos_mode)^2 + sin_mode^2))  # real part
            eig_weak_poisson[n_dofs-k+1] = 0.0  # imaginary part

        end

        # N/2 mode
        coef0 = mass_0[1]
        coef1 = mass_1[1]
        for j=1:s_deg_0 - 1
           coef0 = coef0 + 2 * mass_0[j+1]*cos(pi*j)
           coef1 = coef1 + 2 * mass_1[j+1]*cos(pi*j)
        end
        
        # add last term for larger matrix
        j = s_deg_0
        coef0 = coef0 + 2 * mass_0[j+1]*cos(pi*j)

        # compute eigenvalues
        eig_mass0[n_dofs÷2+1] = coef0
        eig_mass1[n_dofs÷2+1] = coef1
        eig_weak_ampere[n_dofs÷2+1] = 2.0 * (coef1 / coef0)
        eig_weak_poisson[n_dofs÷2+1] = 1.0 / (coef1 * 4.0) 

        new( Lx, delta_x, n_dofs, s_deg_0, s_deg_1, mass_0, mass_1,
             eig_mass0, eig_mass1, eig_weak_ampere, eig_weak_poisson,
             plan_fw, plan_bw, work, wsave )


    end

end 

export compute_rhs_from_function!

"""
   compute_rhs_from_function(self, func, degree, coefs_dofs)

Compute the FEM right-hand-side for a given function f and periodic splines of given degree.

Its components are ``\\int f N_i dx`` where ``N_i`` is the B-spline starting at ``x_i``. 
"""
function compute_rhs_from_function!( coefs_dofs :: Vector{Float64},
                                     self       :: Maxwell1DFEM, 
                                     func       :: Function, 
                                     degree     :: Int64 )

    bspl = zeros(Float64, (degree+1,degree+1))

    # take enough Gauss points so that projection is exact for splines of degree deg
    # rescale on [0,1] for compatibility with B-splines
    x, w = gausslegendre(degree+1)

    x .= 0.5 .* (x .+ 1.0)
    w .= 0.5 .* w

    # Compute bsplines at gauss_points
    for k=1:degree+1
        bspl[k,:] .=  uniform_bsplines_eval_basis(degree, x[k])
    end

    # Compute coefs_dofs = int f(x)N_i(x) 
    for i = 1:self.n_dofs
        coef = 0.0
        # loop over support of B spline
        for j = 1:degree+1
           # loop over Gauss points
            for k=1:degree+1
                coef = coef + w[k]*func(self.delta_x*(x[k] + i + j - 2)) * bspl[k,degree+2-j]
            end
        end
        # rescale by cell size
        coefs_dofs[i] = coef * self.delta_x
     end

end 

function solve_circulant!(self, eigvals, rhs)

    n = self.n_dofs
    # Compute res from rhs, using eigenvalue of circulant  matrix
    self.work .= rhs
    # Forward FFT
    mul!(self.wsave, self.plan_fw, self.work )
    self.wsave[1] = self.wsave[1] * eigvals[1]
    for k=2:n÷2
       re_p = self.wsave[k] * eigvals[k] - self.wsave[n-k+2] * eigvals[n-k+2]
       im_p = self.wsave[k] * eigvals[n-k+2] + self.wsave[n-k+2] * eigvals[k]
       self.wsave[k]     = re_p
       self.wsave[n-k+2] = im_p
    end
    self.wsave[n÷2+1] = self.wsave[n÷2+1]*eigvals[n÷2+1]
    # Backward FFT 
    mul!( self.work, self.plan_bw, self.wsave)

    self.work ./= n

end

export compute_e_from_rho!

function compute_e_from_rho!(e    :: Vector{Float64}, 
                             self :: Maxwell1DFEM, 
                             rho  :: Vector{Float64} )

    # Compute potential phi from rho, using eigenvalue of circulant 
    # inverse matrix
    solve_circulant!(self, self.eig_weak_poisson, rho)
    # Compute spline coefficients of Ex from those of phi
    for i=2:self.n_dofs
        e[i] = self.work[i-1] -  self.work[i]
    end
    # treat Periodic point
    e[1] = self.work[self.n_dofs] - self.work[1] 

end

export compute_e_from_j!
"""
    compute_e_from_j!(e, maxwell_solver, current, component)

Compute ``E_i`` from ``j_i`` integrated over the time interval using weak Ampere formulation
"""
function compute_e_from_j!(e         :: Vector{Float64}, 
                           self      :: Maxwell1DFEM, 
                           current   :: Vector{Float64}, 
                           component :: Int64)

     n = self.n_dofs
     eigvals = zeros(Float64, n)

     # Multiply by inverse mass matrix  using the eigenvalues of the circulant 
     # inverse matrix

     if (component == 1)
         for i=1:n÷2+1
            eigvals[i] = 1.0 / self.eig_mass1[i]
         end
         solve_circulant!(self, eigvals, current)
     elseif (component == 2)
         for i=1:n÷2+1
            eigvals[i] = 1.0 / self.eig_mass0[i]
         end
         solve_circulant!(self, eigvals, current)
     else
         throw(ArgumentError("Component $component not implemented "))
     end
     
     # Update the electric field and scale
     e .= e - self.work ./ self.delta_x

end

export l2norm_squared

"""
    l2norm_squared(maxwell_solver, coefs_dofs, degree)

Compute square of the L2norm 

"""
function l2norm_squared(self, coefs_dofs, degree)

    # Multiply coefficients by mass matrix (use diagonalization FFT and mass matrix eigenvalues)
    if (degree == self.s_deg_0 )

        solve_circulant!(self, self.eig_mass0, coefs_dofs)

    elseif (degree == self.s_deg_1)

        solve_circulant!(self, self.eig_mass1, coefs_dofs)

    end

    # Multiply by the coefficients from the left (inner product)
    r = sum(coefs_dofs .* self.work)
    # Scale by delt_x
    r .* self.delta_x

end 



export l2norm_squared2

"""
    l2norm_squared(maxwell_solver, coefs_dofs, degree)

Compute square of the L2norm 

"""
function l2norm_squared2(self, coefs_dofs, degree)

    # Multiply coefficients by mass matrix (use diagonalization FFT and mass matrix eigenvalues)
    if (degree == self.s_deg_0 )

        solve_circulant!(self, self.eig_mass0, coefs_dofs)

    elseif (degree == self.s_deg_1)

        solve_circulant!(self, self.eig_mass1, coefs_dofs)

    end

    # Multiply by the coefficients from the left (inner product)
    r = sum(self.work .* self.work)
    # Scale by delt_x
    r .* self.delta_x

end 




export l2projection!

"""
    l2projection!(coefs_dofs, maxwell, func, degree)

Compute the L2 projection of a given function f on periodic splines 
of given degree
"""
function l2projection!(coefs_dofs :: Vector{Float64},
                       self       :: Maxwell1DFEM, 
                       func       :: Function, 
                       degree     :: Int64)

    n = self.n_dofs
    eigvals = zeros(Float64, n)

    # Compute right-hand-side
    compute_rhs_from_function!( coefs_dofs, self, func, degree)

    # Multiply by inverse mass matrix 
    if (degree == self.s_deg_0)
       for i=1:n÷2+1
          eigvals[i] = 1.0 / self.eig_mass0[i]
       end
    elseif  (degree == self.s_deg_0-1)
       for i=1:n÷2+1
          eigvals[i] = 1.0 / self.eig_mass1[i]
       end
    else
       throw(ArgumentError("degree $degree not available")) 
    end

    solve_circulant!(self, eigvals, coefs_dofs)

    # Account for scaling in the mass matrix by dx

    coefs_dofs .= coefs_dofs ./ self.delta_x

end
  
export compute_e_from_b!

"""
    compute_e_from_b!(field_out, maxwell_solver, delta_t, field_in)

compute Ey from Bz using weak Ampere formulation 

"""
function compute_e_from_b!(field_out :: Vector{Float64}, 
                           self      :: Maxwell1DFEM, 
                           delta_t   :: Float64, 
                           field_in  :: Vector{Float64} )
    
    coef = delta_t / self.delta_x

    # Compute potential weak curl of Bz using eigenvalue of circulant inverse matrix
    solve_circulant!(self, self.eig_weak_ampere, field_in)
    # Update bz from self value
    field_out .+= coef .* self.work

end

export compute_b_from_e!
"""
    compute_b_from_e!( field_out, maxwell_solver, delta_t, field_in) 

Compute Bz from Ey using strong 1D Faraday equation for spline coefficients
```math
B_z^{new}(x_j) = B_z^{old}(x_j) - \\frac{\\Delta t}{\\Delta x} (E_y(x_j) - E_y(x_{j-1})
```
"""
function compute_b_from_e!( field_out :: Vector{Float64},
                            self      :: Maxwell1DFEM, 
                            delta_t   :: Float64, 
                            field_in  :: Vector{Float64}) 

    coef = delta_t/self.delta_x
    # relation betwen spline coefficients for strong Ampere
    for i=2:self.n_dofs
       field_out[i] = field_out[i] + coef * ( field_in[i-1] - field_in[i] )
    end
    # treat Periodic point
    field_out[1] = field_out[1] + coef * ( field_in[end] - field_in[1] )

end




export compute_derivatives_from_basis!

function compute_derivatives_from_basis!( field_out :: Vector{Float64},
                            self      :: Maxwell1DFEM, 
                            field_in  :: Vector{Float64}) 

    coef = 1/self.delta_x
    # relation betwen spline coefficients for strong Ampere
    for i=1:(self.n_dofs-1)
       field_out[i] =  coef * ( field_in[i] - field_in[i+1] )
    end
    # treat Periodic point
    field_out[end] =  coef * ( field_in[end] - field_in[1] )

end

export compute_derivatives_from_basis2!

function compute_derivatives_from_basis2!( field_out :: Vector{Float64},
                            self      :: Maxwell1DFEM, 
                            field_in  :: Vector{Float64}) 

    coef = 1/self.delta_x
    # relation betwen spline coefficients for strong Ampere
    for i=2:(self.n_dofs)
       field_out[i] =  coef * ( field_in[i] - field_in[i-1] )
    end
    # treat Periodic point
    field_out[1] =  coef * ( field_in[1] - field_in[end] )

end



"""
    inner_product( maxwell_solver, coefs1_dofs, coefs2_dofs, degree )

-  `maxwell_solver` : Maxwell solver object
-  `coefs1_dofs` : Coefficient for each DoF
-  `coefs2_dofs` : Coefficient for each DoF
-  `degree : Specify the degree of the basis functions

return squared L2 norm

"""
function inner_product( self, coefs1_dofs, coefs2_dofs, degree ) 

     # Multiply coefficients by mass matrix (use diagonalization FFT and mass matrix eigenvalues)
     if degree == self.s_deg_0

         solve_circulant!(self, self.eig_mass0, coefs2_dofs)

     elseif degree == self.s_deg_1

         solve_circulant!(self, self.eig_mass1, coefs2_dofs)

     end 

     # Multiply by the coefficients from the left (inner product)
     # Scaled by delt_x

     sum(coefs1_dofs .* self.work) * self.delta_x
     
end

