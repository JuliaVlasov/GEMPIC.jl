abstract type AbstractGrid end

export OneDGrid, TwoDGrid

"""
    TwoDGrid( dimx, nx, dimy, ny)

Generate a cartesians mesh on rectangle `dimx`x `dimy` with `nx` x `ny` points

- `nx` : indices are in [1:nx]
- `ny` : indices are in [1:ny]
- `dimx = xmax - xmin`
- `dimy = ymax - ymin`
- `x, y` : node positions
- `dx, dy` : step size
"""
struct TwoDGrid <: AbstractGrid

    nx::Int
    ny::Int
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    dimx::Float64
    dimy::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    dx::Float64
    dy::Float64

    function TwoDGrid(xmin, xmax, nx :: Int, ymin, ymax, ny :: Int)

        dimx = xmax - xmin
        dimy = ymax - ymin

        x = LinRange(0, dimx, nx + 1) |> collect
        y = LinRange(0, dimy, ny + 1) |> collect

        dx = dimx / nx
        dy = dimy / ny

        new(nx, ny, xmin, xmax, ymin, ymax, dimx, dimy, x, y, dx, dy)

    end

end

TwoDGrid(dimx, nx::Int, dimy, ny::Int) = TwoDGrid(0.0, dimx, nx, 0.0, dimy, ny)

"""
    TwoDGrid( xmin, xmax, nx, ymin, ymax, ny )

Simple structure to store mesh data from 2 dimensions
"""
struct OneDGrid <: AbstractGrid

    nx::Int
    xmin::Float64
    xmax::Float64
    dimx::Float64
    dx::Float64
    x::Vector{Float64}

    function OneDGrid(xmin, xmax, nx::Int)

        dimx = xmax - xmin
        dx = dimx / (nx - 1)
        x = LinRange(0, dimx, nx + 1) |> collect

        new(nx, xmin, xmax, dimx, dx, x)

    end

end

export get_x

"""  
    get_x( mesh, i )

Get position
"""
function get_x(m::OneDGrid, i)

    m.xmin + (i - 1) * m.dx

end

"""  
    get_cell_and_offset( mesh, x )

Get cell and offset

We compute the cell indices where the particle is and its relative 
normalized position inside the cell

"""
function get_cell_and_offset(m::OneDGrid, x)

    cell = floor(Int64, ((x - m.xmin) / m.Lx) * m.nx) + 1
    offset = (x - get_x(m, cell)) / dx

    cell, offset

end
