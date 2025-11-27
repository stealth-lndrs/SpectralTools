"""
    cheb(N)

Return the Chebyshev–Lobatto points `x` and the differentiation matrix `D`
for order `N`, following Trefethen's implementation.
"""
function cheb(N::Integer)
    N < 0 && throw(ArgumentError("N must be non-negative"))
    if N == 0
        return ones(1), zeros(1, 1)
    end

    x = cos.(pi .* (0:N) ./ N)
    c = [2.0; ones(N - 1); 2.0] .* (-1).^(0:N)
    xcol = reshape(x, N + 1, 1)
    X = repeat(xcol, 1, N + 1)
    dX = X .- X'
    D = (c * (1 ./ c)') ./ (dX + I)
    D .-= Diagonal(sum(D, dims = 2)[:])
    return x, D
end

"""
    cheb_lobatto_weights(N)

Return Chebyshev–Lobatto quadrature weights as described in the slides.
"""
function cheb_lobatto_weights(N::Integer)
    N <= 0 && throw(ArgumentError("N must be positive"))
    weights = zeros(Float64, N + 1)
    invN2 = 1.0 / (N^2)
    weights[1] = invN2
    weights[end] = invN2
    for j in 2:N
        if isodd(j - 1)
            weights[j] = 2 * invN2
        else
            weights[j] = 0.0
        end
    end
    return weights
end
