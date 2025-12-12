#!/usr/bin/env julia
#
# Usage Example 4: Nonlinear BVP y'' = sin(y) with homogeneous Dirichlet BCs.
# Demonstrates solve_nonlinear_bvp with analytic derivatives and reports
# Newton convergence diagnostics.

using PoliSpectralTools
using LinearAlgebra: norm
using DelimitedFiles: writedlm
using Printf

g(x, y, yp) = sin(y)
dgdy(x, y, yp) = cos(y)
dg_dyp(x, y, yp) = zero(x)

result = solve_nonlinear_bvp(g;
    dg_dy = dgdy,
    dg_dyp = dg_dyp,
    N = 48,
    basis = :chebyshev,
    domain = (-1.0, 1.0),
    bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.7)),
    initial_guess = x -> 0.35 * (x + 1),  # linear ramp matching boundary data
    tol = 1e-10,
    maxiter = 20,
)

grid = result.grid
residual = grid.D2 * result.u .- g.(grid.x, result.u, grid.D1 * result.u)
@printf("Nonlinear BVP iterations = %d, residual norm = %.3e (converged = %s)\n",
        result.iterations, norm(residual, Inf), result.converged)

out_path = joinpath(@__DIR__, "..", "docs", "assets", "bvp_nonlinear_solution.csv")
writedlm(out_path, [result.x result.u])
