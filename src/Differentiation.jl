"""
    chebD2(N)

Return the interior block of the squared Chebyshev differentiation matrix
alongside the interior grid points.
"""
function chebD2(N::Integer)
    N < 2 && throw(ArgumentError("N must be at least 2"))
    x, D = cheb(N)
    D2 = D * D
    inner = 2:N
    return D2[inner, inner], x[inner]
end
