export Mesh

"""
    Mesh( xmin, xmax,nx )

Simple structure to store mesh data from 1 to 3 dimensions
"""
struct Mesh 

    nx    :: Vector{Int}
    xmin  :: Vector{Real}
    xmax  :: Vector{Real}
    Lx    :: Vector{Real} 
    dx    :: Vector{Real}

    function Mesh( nx, xmin, xmax, Lx, dx)

        new( nx, xmin, xmax, Lx, dx )

    end

    function Mesh( xmin :: Real, xmax :: Real, nx :: Int) 

        Lx = xmax - xmin
        dx = Lx / (nx - 1)

        Mesh( [nx], [xmin], [xmax], [Lx], [dx] )

    end

    function Mesh( x1min :: Real, x1max :: Real, nx1 :: Int,
                   x2min :: Real, x2max :: Real, nx2 :: Int  )

        nx   = [nx1, nx2]
        xmin = [x1min, x2min]
        xmax = [x1max, x2max]
        Lx   = xmax .- xmin
        dx   = Lx ./ ( nx .- 1 )

        Mesh( nx, xmin, xmax, Lx, dx )

    end

    function Mesh( x1min :: Real, x1max :: Real, nx1 :: Int,
                   x2min :: Real, x2max :: Real, nx2 :: Int,
                   x3min :: Real, x3max :: Real, nx3 :: Int  )

        nx   = [nx1, nx2, nx3]
        xmin = [x1min, x2min, x3min]
        xmax = [x1max, x2max, x3max]
        Lx   = xmax .- xmin
        dx   = Lx ./ ( nx .- 1 )

        Mesh( nx, xmin, xmax, Lx, dx )

    end

end 

export get_x

"""  
    get_x( mesh, i )

Get position
"""
function get_x( m :: Mesh, i )

    m.xmin .+ (i .- 1) .* m.dx
    
end 

"""  
    get_cell_and_offset( mesh, x )

Get cell and offset

We compute the cell indices where the particle is and its relative 
normalized position inside the cell

"""
function get_cell_and_offset( m :: Mesh, x ) :: Int64

    cell   = floor.(Int64, ((x .- m.xmin) ./ m.Lx) .* m.nx) .+ 1
    offset = (x .- get_x( m, cell)) ./ dx

	cell, offset
    
end 
