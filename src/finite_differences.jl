export
    CFDM_3_1, CFDM_3_2,
    finitediff!, laplacian!, divergence!

# pre-computed finite difference stencils
global const CFDM_3_1 = central_fdm(3, 1)
global const CFDM_3_2 = central_fdm(3, 2)

# derivative functions for internal usage
function ∂x(u::AbstractVector{Float64}, i::Int, FDM)
    s::Float64 = 0.0
    for (q,h) in zip(FDM.grid, FDM.coefs)
        s += u[i+q] * h
    end # for
    return s
end # function

function ∂x(u::AbstractMatrix{Float64}, i::Int, j::Int, FDM)
    s::Float64 = 0.0
    for (q,h) in zip(FDM.grid, FDM.coefs)
        s += u[i+q, j] * h
    end # for
    return s
end # function

function ∂y(u::AbstractMatrix{Float64}, i::Int, j::Int, FDM)
    s::Float64 = 0.0
    for (q,h) in zip(FDM.grid, FDM.coefs)
        s += u[i, j+q] * h
    end # for
    return s
end # function

function ∂x(u::AbstractArray{Float64,3}, i::Int, j::Int, k::Int, FDM)
    s::Float64 = 0.0
    for (q,h) in zip(FDM.grid, FDM.coefs)
        s += u[i+q, j, k] * h
    end # for
    return s
end # function

function ∂y(u::AbstractArray{Float64,3}, i::Int, j::Int, k::Int, FDM)
    s::Float64 = 0.0
    for (q,h) in zip(FDM.grid, FDM.coefs)
        s += u[i, j+q, k] * h
    end # for
    return s
end # function

function ∂z(u::AbstractArray{Float64,3}, i::Int, j::Int, k::Int, FDM)
    s::Float64 = 0.0
    for (q,h) in zip(FDM.grid, FDM.coefs)
        s += u[i, j, k+q] * h
    end # for
    return s
end # function


"""
    finitediff!(ux, u, a, FDM=CFDM_3_1)
    finitediff!(ux, uy, u, ax, ay, FDM=CFDM_3_1)
    finitediff!(ux, uy, u, a, FDM=CFDM_3_1)
    finitediff!(ux, uy, uz, u, ax, ay, az, FDM=CFDM_3_1)
    finitediff!(ux, uy, uz, u, a, FDM=CFDM_3_1)
Evaluate derivative of `u` through the finite difference stencil `FDM`
and store the results in the appropriate arrays.

If `u` has more than one dimension, it is assumed that the dimensions are
ordered as (x, y, z).
It is required that `size(ux) == size(uy) == size(uz) == size(u)`.

By default, `FDM` is set to the 3-point stencil of the first derivative.
`ax`, `ay` and `az` are the appropriate discretization factors
(1/dx^n, 1/dy^n, 1/dx^n for the n-th derivative);
if a single value `a` is given, it is expanded so that
`ax = ay = az = a`.
"""
@views function finitediff!(ux::A, u::A, a::Real, FDM=CFDM_3_1) where
    {A<:AbstractVector{Float64}}
    n = length(u)
    i_start = 1 + abs(FDM.grid[1])
    i_end = n - FDM.grid[end]
    for i in i_start:i_end
        ux[i] = a * ∂x(u,i,FDM)
    end # for
end # function

@views function finitediff!(ux::A, uy::A, u::A,
                            ax::B, ay::B, FDM=CFDM_3_1) where
    {A<:AbstractMatrix{Float64}, B<:Real}
    nx, ny = size(u)
    i_start = j_start = 1 + abs(FDM.grid[1])
    i_end = nx - FDM.grid[end]
    j_end = ny - FDM.grid[end]
    for j in j_start:j_end, i in i_start:i_end
        ux[i,j] = ax * ∂x(u,i,j,FDM)
        uy[i,j] = ay * ∂y(u,i,j,FDM)
    end # for
end # function

finitediff!(ux::A, uy::A, u::A, a::Real, FDM=CFDM_3_1) where
    {A<:AbstractMatrix{Float64}} = finitediff!(ux, uy, u, a, a, FDM)

@views function finitediff!(ux::A, uy::A, uz::A, u::A,
                            ax::B, ay::B, az::B, FDM=CFDM_3_1) where
    {A<:AbstractArray{Float64,3}, B<:Real}
    nx, ny, nz = size(u)
    i_start = j_start = k_start = 1 + abs(FDM.grid[1])
    i_end = nx - FDM.grid[end]
    j_end = ny - FDM.grid[end]
    k_end = nz - FDM.grid[end]
    for k in k_start:k_end, j in j_start:j_end, i in i_start:i_end
        ux[i,j,k] = ax * ∂x(u,i,j,k,FDM)
        uy[i,j,k] = ay * ∂y(u,i,j,k,FDM)
        uz[i,j,k] = az * ∂z(u,i,j,k,FDM)
    end # for
end # function

finitediff!(ux::A, uy::A, uz::A, u::A, a::Real, FDM=CFDM_3_1) where
    {A<:AbstractArray{Float64,3}} = finitediff!(ux, uy, uz, u, a, a, a, FDM)


"""
    laplacian!(du, u, a, FDM=CFDM_3_2)
    laplacian!(du, u, ax, ay, FDM=CFDM_3_2)
    laplacian!(du, u, ax, ay, az, FDM=CFDM_3_2)
Evaluate the finite-difference laplacian of `u` (∇²u),
storing the result in `du`.
By default, the 3-point stencil of the 2nd derivative is used.
"""
@views function laplacian!(du::A, u::A, a::B, FDM=CFDM_3_2) where
    {A<:AbstractVector{Float64}, B<:Real}
    finitediff!(du, u, a, FDM)
end # function

@views function laplacian!(du::A, u::A, ax::B, ay::B, FDM=CFDM_3_2) where
    {A<:AbstractMatrix{Float64}, B<:Real}
    nx, ny = size(u)
    i_start = j_start = 1 + abs(FDM.grid[1])
    i_end = nx - FDM.grid[end]
    j_end = ny - FDM.grid[end]
    for j in j_start:j_end, i in i_start:i_end
        du[i,j] = ax*∂x(u,i,j,FDM) + ay*∂y(u,i,j,FDM)
    end # for
end # function

@views function laplacian!(du::A, u::A, ax::B, ay::B, az::B, FDM=CFDM_3_2) where
    {A<:AbstractArray{Float64,3}, B<:Real}
    nx, ny, nz = size(u)
    i_start = j_start = k_start = 1 + abs(FDM.grid[1])
    i_end = nx - FDM.grid[end]
    j_end = ny - FDM.grid[end]
    k_end = nz - FDM.grid[end]
    for k in k_start:k_end, j in j_start:j_end, i in i_start:i_end
        du[i,j,k] = ax*∂x(u,i,j,k,FDM) + ay*∂y(u,i,j,k,FDM) + az*∂z(u,i,j,k,FDM)
    end # for
end

laplacian!(du, u, a::Real, FDM=CFDM_3_2) =
    laplacian!(du, u, ntuple(_->a, ndims(u))..., FDM)

"""
    divergence!(du, u, a, FDM=CFDM_3_1)
    divergence!(du, ux, uy, ax, ay, FDM=CFDM_3_1)
    divergence!(du, ux, uy, uz, ax, ay, az, FDM=CFDM_3_1)
Evaluate the finite-difference divergence of `u` (∇⋅u),
storing the result in `du`.
By default, the 3-point stencil of the 1st derivative is used.
"""
divergence!(du::A, u::A, a::Real, FDM=CFDM_3_1) where
    {A<:AbstractVector{Float64}} = laplacian!(du, u, a, FDM)

@views function divergence!(du::A, ux::A, uy::A,
                            ax::B, ay::B, FDM=CFDM_3_1) where
    {A<:AbstractMatrix{Float64}, B<:Real}
    nx, ny = size(du)
    i_start = j_start = 1 + abs(FDM.grid[1])
    i_end = nx - FDM.grid[end]
    j_end = ny - FDM.grid[end]
    for j in j_start:j_end, i in i_start:i_end
        du[i,j] = ax*∂x(ux,i,j,FDM) + ay*∂y(uy,i,j,FDM)
    end # for
end # function

@views function divergence!(du::A, ux::A, uy::A, uz::A,
                            ax::B, ay::B, az::B, FDM=CFDM_3_1) where
    {A<:AbstractArray{Float64,3}, B<:Real}
    nx, ny, nz = size(du)
    i_start = j_start = k_start = 1 + abs(FDM.grid[1])
    i_end = nx - FDM.grid[end]
    j_end = ny - FDM.grid[end]
    k_end = nz - FDM.grid[end]
    for k in k_start:k_end, j in j_start:j_end, i in i_start:i_end
        du[i,j,k] = ax*∂x(ux,i,j,k,FDM) + ay*∂y(uy,i,j,k,FDM) + az*∂z(uz,i,j,k,FDM)
    end # for
end
