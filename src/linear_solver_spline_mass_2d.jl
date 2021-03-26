struct TwoDLinearSolverSplineMass

    nx::Int
    ny::Int
    factor::Float64
    eig_values_1::Vector{Float64}
    eig_values_2::Vector{Float64}

    function TwoDLinearSolverSplineMass(nx, ny, eig_values_1, eig_values_2)

        factor = 1.0
        new(nx, ny, factor, eig_values_1, eig_values_2)

    end

end

function solve_real_mass1(solver::TwoDLinearSolverSplineMass, rhs)

    nx, ny = solver.nx, solver.ny
    array1d_x = zeros(ComplexF64, nx)
    array1d_y = zeros(ComplexF64, ny)
    scratch = zeros(ComplexF64, nx, ny)
    k = 0
    for j = 1:ny
        for i = 1:nx
            k = k + 1
            array1d_x[i] = rhs[k]
        end
        fft!(array1d_x)
        scratch[:, j] .= array1d_x
    end
    for i = 1:nx
        array1d_y .= scratch[i, :]
        fft!(array1d_y)
        scratch[i, :] .= array1d_y
    end

    # Multiply by inverse mass
    for j = 1:ny, i = 1:nx
        scratch[i, j] /= (solver.eig_values_1[i] * solver.eig_values_2[j] * solver.factor)
    end

    for i = 1:nx
        array1d_y .= scratch[i, :]
        fft!(array1d_y)
        scratch[i, :] .= array1d_y
    end

    sol = zero(rhs)
    k = 0
    for j = 1:ny
        array1d_x .= scratch[:, j]
        ifft!(array1d_x)
        for i = 1:nx
            k = k + 1
            sol[k] = real(array1d_x[i])
        end
    end

end
