"""
    poisson_chebyshev_2d(f_fun, N)

Solve the Poisson problem Δu = f on [-1, 1]^2 with homogeneous boundary
conditions using Chebyshev–Lobatto collocation of order `N`.
"""
function poisson_chebyshev_2d(f_fun, N::Integer)
    N < 2 && throw(ArgumentError("N must be at least 2"))
    x, D = cheb(N)
    D2 = D * D
    inner = 2:N
    D2_inner = D2[inner, inner]
    x_inner = x[inner]
    nint = length(x_inner)
    In = Matrix{Float64}(I, nint, nint)
    L = kron(In, D2_inner) + kron(D2_inner, In)

    xx = repeat(reshape(x_inner, 1, nint), nint, 1)
    yy = repeat(reshape(x_inner, nint, 1), 1, nint)
    fvals = f_fun.(xx, yy)
    rhs = vec(fvals)
    u_inner = L \ rhs

    U = zeros(nint + 2, nint + 2)
    U[inner, inner] .= unvec(u_inner, nint, nint)
    return U, x, x
end
