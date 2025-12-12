module Quadrature

using ..Chebyshev: ChebyQuadrature

export chebyshev_quadrature, chebyshev_gauss_integral, chebyshev_lobatto_integral

"""
    chebyshev_quadrature(f, N; domain = (-1.0, 1.0), kind = :gauss)

Evaluate the integral of `f` over `domain` by Chebyshev quadrature using `N`
nodes.  `kind` accepts `:gauss`, `:lobatto`, `:radau_left`, or `:radau_right`.
Returns a tuple `(integral, nodes, weights)`.
"""
function chebyshev_quadrature(f::Function, N::Integer;
        domain::Tuple = (-1.0, 1.0), kind::Symbol = :gauss)
    tipo = _kind_to_tipo(kind)
    nodes, weights = ChebyQuadrature(domain[1], domain[2], N, tipo)
    vals = f.(nodes)
    return (sum(weights .* vals), nodes, weights)
end

"""
    chebyshev_gauss_integral(f, N; domain = (-1.0, 1.0))

Convenience wrapper around `chebyshev_quadrature` with Gauss nodes.
"""
chebyshev_gauss_integral(f::Function, N::Integer; domain::Tuple = (-1.0, 1.0)) =
    first(chebyshev_quadrature(f, N; domain = domain, kind = :gauss))

"""
    chebyshev_lobatto_integral(f, N; domain = (-1.0, 1.0))

Convenience wrapper using Lobatto (Clenshaw–Curtis) nodes.
"""
chebyshev_lobatto_integral(f::Function, N::Integer; domain::Tuple = (-1.0, 1.0)) =
    first(chebyshev_quadrature(f, N; domain = domain, kind = :lobatto))

function _kind_to_tipo(kind::Symbol)
    kind === :gauss && return 1
    kind === :lobatto && return 2
    kind === :radau_left && return 3
    kind === :radau_right && return 4
    error("Unknown quadrature kind $kind. Use :gauss, :lobatto, :radau_left, or :radau_right.")
end

end
