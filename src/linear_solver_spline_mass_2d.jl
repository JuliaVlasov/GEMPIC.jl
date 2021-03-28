struct TwoDLinearSolverSplineMass

    nx1::Int
    nx2::Int
    factor::Float64
    eigvals1::Vector{Float64}
    eigvals2::Vector{Float64}
    wk::Array{ComplexF64,2}

    function TwoDLinearSolverSplineMass(nx1, nx2, eigvals1, eigvals2)

        factor = 1.0
        wk = zeros(ComplexF64, nx1, nx2)
        new(nx1, nx2, factor, eigvals1, eigvals2, wk)

    end

end

function solve(solver::TwoDLinearSolverSplineMass, rhs::Vector{Float64})::Vector{Float64}

    nx1, nx2 = solver.nx1, solver.nx2

    solver.wk .= reshape(rhs, nx1, nx2)

    fft!(solver.wk)

    for j = 1:nx2, i = 1:nx1
        @inbounds solver.wk[i, j] /=
            (solver.eigvals1[i] * solver.eigvals2[j] * solver.factor)
    end

    ifft!(solver.wk)

    vec(real(solver.wk))

end
