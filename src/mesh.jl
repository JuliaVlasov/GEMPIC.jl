abstract type AbstractGrid end

export OneDGrid, TwoDGrid, ThreeDGrid

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

    function TwoDGrid(xmin, xmax, nx::Int, ymin, ymax, ny::Int)
        dimx = xmax - xmin
        dimy = ymax - ymin

        x = collect(LinRange(xmin, xmax, nx + 1))
        y = collect(LinRange(ymin, ymax, ny + 1))

        dx = dimx / nx
        dy = dimy / ny

        return new(nx, ny, xmin, xmax, ymin, ymax, dimx, dimy, x, y, dx, dy)
    end
end

TwoDGrid(dimx, nx::Int, dimy, ny::Int) = TwoDGrid(0.0, dimx, nx, 0.0, dimy, ny)

"""
    OneDGrid( xmin, xmax, nx )

Simple structure to store mesh data from 1 dimension
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
        x = collect(LinRange(xmin, xmax, nx + 1))

        return new(nx, xmin, xmax, dimx, dx, x)
    end
end

export get_x

"""  
    get_x( mesh, i )

Get position
"""
function get_x(m::OneDGrid, i)
    return m.xmin + (i - 1) * m.dx
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

    return cell, offset
end

"""
    ThreeDGrid( dimx, nx, dimy, ny, dimz, nz)

Generate a cartesians mesh on cube `dimx` x `dimy` x `dimz` with `nx` x `ny` x `nz` points

- `nx` : indices are in [1:nx]
- `ny` : indices are in [1:ny]
- `nz` : indices are in [1:nz]
- `dimx = xmax - xmin`
- `dimy = ymax - ymin`
- `dimz = zmax - zmin`
- `x, y, z` : node positions
- `dx, dy, dz` : step size
"""
struct ThreeDGrid <: AbstractGrid
    nx::Int
    ny::Int
    nz::Int
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    zmin::Float64
    zmax::Float64
    dimx::Float64
    dimy::Float64
    dimz::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    dx::Float64
    dy::Float64
    dz::Float64

    function ThreDGrid(xmin, xmax, nx::Int, ymin, ymax, ny::Int, zmin, zmax, nz::Int)
        dimx = xmax - xmin
        dimy = ymax - ymin
        dimz = zmax - zmin

        x = collect(LinRange(xmin, xmax, nx + 1))
        y = collect(LinRange(ymin, ymax, ny + 1))
        z = collect(LinRange(zmin, zmax, nz + 1))

        dx = dimx / nx
        dy = dimy / ny
        dz = dimz / nz

        return new(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, dimx, dimy, dimz, x, y, z, dx, dy, dz)
    end
end

ThreeDGrid(dimx, nx::Int, dimy, ny::Int, dimz, nz::Int) = ThreDGrid(0.0, dimx, nx, 0.0, dimy, ny, 0.0, dimz, nz)
