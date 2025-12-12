#!/usr/bin/env julia
#
# Quadrature example: approximate definite integrals via Chebyshev rules.

using PoliSpectralTools
using Printf

f(x) = exp(x)
domain = (0.0, 1.0)
exact = exp(1) - 1

I_gauss, nodes_g, weights_g = chebyshev_quadrature(f, 32; domain = domain, kind = :gauss)
@printf("Gauss-Chebyshev integral exp(x) dx on [0,1]  -> %.10f (error %.2e)\n",
        I_gauss, abs(I_gauss - exact))

g(x) = cos(3x)
domain2 = (-1.0, 1.0)
exact2 = (sin(3) - sin(-3)) / 3

I_lobatto = chebyshev_lobatto_integral(g, 48; domain = domain2)
@printf("Lobatto-Chebyshev integral cos(3x) dx on [-1,1] -> %.10f (error %.2e)\n",
        I_lobatto, abs(I_lobatto - exact2))
