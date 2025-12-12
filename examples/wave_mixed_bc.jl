#!/usr/bin/env julia
#
# Usage Example 5: Wave equation with mixed Neumann/Dirichlet boundary conditions.

using PoliSpectralTools
using PoliSpectralTools.Collocation: build_grid
using LinearAlgebra: norm
using Printf
using Plots: plot, @animate, gif

c = 1.0
u0(x) = sin(π * (x + 1) / 2)
v0(x) = zero(x)
flux(t) = cos(5t)

bc = (left = (:neumann, (x, t) -> flux(t)), right = (:dirichlet, 0.0))

sol = solve_wave_1d(u0, v0, (0.0, 2.0);
    N = 50,
    dt = 1e-3,
    c = c,
    basis = :chebyshev,
    domain = (-1.0, 1.0),
    bc = bc,
)

grid = build_grid(50; basis = :chebyshev, domain = (-1.0, 1.0))
function energy(u_slice, v_half)
    grad = grid.D1 * u_slice
    return 0.5 * (norm(v_half)^2 + c^2 * norm(grad)^2)
end

E0 = energy(sol.u[:, 1], sol.v[:, 1])
Eend = energy(sol.u[:, end], sol.v[:, end])
drift = abs(Eend - E0) / E0

left_flux = (grid.D1[1, :] * sol.u[:, end])
@printf("Wave mixed-BC run finished. Energy drift = %.3e, left flux error = %.3e\n",
        drift, left_flux - flux(sol.t[end]))

# Generate GIF showing reflections
ENV["GKSwstype"] = "100"
anim = @animate for k in 1:length(sol.t)
    plot(sol.x, sol.u[:, k];
         ylim = (-0.8, 0.8),
         xlabel = "x",
         ylabel = "u(x,t)",
         title = @sprintf("Mixed BC standing wave, t = %.2f", sol.t[k]),
         legend = false,
         color = :cyan,
         background_color = :black,
         foreground_color = :white,
         lw = 2)
end
gif(anim, joinpath(@__DIR__, "..", "docs", "assets", "wave_reflection.gif"), fps = 20)
