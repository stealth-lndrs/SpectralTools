#!/usr/bin/env julia
#
# Demonstrates solve_nonlinear_bvp (Newton) on y'' + y^3 = cos(πx)

using PoliSpectralTools
using LinearAlgebra: norm
using Printf

g(x, y, yp) = cospi(x) - y^3
dg_dy(x, y, yp) = -3y^2
dg_dyp(x, y, yp) = zero(x)

sol = solve_nonlinear_bvp(g;
    dg_dy = dg_dy,
    dg_dyp = dg_dyp,
    N = 60,
    basis = :chebyshev,
    domain = (-1.0, 1.0),
    bc = (left = (:dirichlet, 0.2), right = (:dirichlet, -0.1)),
    initial_guess = x -> 0.2 - 0.15 * (x + 1),
    tol = 1e-10,
    maxiter = 25,
)

residual = sol.grid.D2 * sol.u .- g.(sol.x, sol.u, sol.grid.D1 * sol.u)
@printf("Newton solve (y'' + y^3 = cos(πx)) iterations = %d, residual = %.2e\n",
        sol.iterations, norm(residual, Inf))
