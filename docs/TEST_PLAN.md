---
layout: default
title: "Test Plan & Example Roadmap"
---

# Test Plan & Example Roadmap

This page summarizes the outstanding work described in
`EXAMPLES_TEST_PLAN.md` (root of the repository) and mirrors the same guidance
inline so contributors do not need to leave the GitHub Pages site.  Pair this
with [Implementation Status](IMPLEMENTATION_STATUS.md) for more detailed steps.

## Usage Examples

| # | Title | Status | Summary |
| --- | --- | --- | --- |
| 1 | Chebyshev linear BVP | ✅ | Implemented (`examples/bvp_linear.jl`) and documented with residual analysis. |
| 2 | Diffusion decay | ✅ | Implemented (`examples/diffusion.jl`) with analytic comparison + GIF. |
| 3 | Legendre vs Chebyshev | ✅ | Implemented (`examples/bvp_legendre.jl`) with error table. |
| 4 | Nonlinear BVP (`y'' = sin y`) | ✅ | Implemented (`examples/bvp_nonlinear.jl`). |
| 5 | Wave (Neumann/Dirichlet) | ✅ | Implemented (`examples/wave_mixed_bc.jl`) and animated. |
| 6 | Wave traveling pulse | ⬜ | **Missing:** Gaussian pulse example with CFL study and optional sponge. See slides PE_Aula_10_N.pdf. |
| 7 | Forced diffusion | ⬜ | **Missing:** Source-driven heat equation tracking injected energy. Slides PE_Aula_06_N.pdf. |
| 8 | Poisson on square | ⬜ | **Missing:** Manufactured solution with convergence table (PE_Aula_09_N.pdf). |
| 9 | Poisson on rectangle | ⬜ | **Missing:** Anisotropic domain showing scaling in `Dx²`/`Dy²`. |
| 10 | Mapping prototype (optional) | ⬜ | **Missing:** Boundary-layer mapping once `Mapping.jl` utility exists. |

## Julia Tests

| # | Description | Status | Notes |
| --- | --- | --- | --- |
| 1 | Chebyshev BVP residual | ✅ | In `Usage-driven Tests`. |
| 2 | Legendre grid consistency | ✅ | Included already. |
| 3 | Nonlinear BVP convergence | ✅ | Included already. |
| 4 | Diffusion analytic comparison | ✅ | Included already. |
| 5 | Wave energy check | ✅ | Included already. |
| 6 | Wave pulse CFL study | ⬜ | Add stable vs unstable runs once Example 6 exists. |
| 7 | Diffusion forcing balance | ⬜ | Validate integral tracking for Example 7. |
| 8 | Poisson square accuracy | ⬜ | Manufactured solution error check for Example 8. |
| 9 | Poisson rectangular scaling | ⬜ | Verify scaling and error target for Example 9. |
| 10 | Mapping impact | ⬜ | Compare mapped vs unmapped errors when mapping utilities arrive. |

## How to Contribute

1. Pick an item marked ⬜ (coordinate via issues/PRs).
2. Follow the detailed checklist in [Implementation Status](IMPLEMENTATION_STATUS.md).
3. Update `EXAMPLES_TEST_PLAN.md` and this summary if statuses change.
4. Run `test/runtests.jl` and regenerate the [test report](TEST_REPORTS.md)
   before opening a pull request.
