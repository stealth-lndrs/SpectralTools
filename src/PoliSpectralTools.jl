
module PoliSpectralTools

include("Generic.jl")
include("Chebyshev.jl")
include("Legendre.jl")
include("Fourier.jl")
include("BoundaryConditions.jl")
include("Collocation.jl")
include("Quadrature.jl")
include("BVP.jl")
include("PDE.jl")

using .Generic
using .Chebyshev
using .Legendre
using .Fourier
using .BoundaryConditions
using .Collocation
using .Quadrature
using .BVP
using .PDE

export Generic, Chebyshev, Legendre, Fourier,
       BoundaryConditions, Collocation, Quadrature,
       solve_linear_bvp, solve_nonlinear_bvp,
        solve_diffusion_1d, solve_wave_1d, solve_poisson_2d,
        chebyshev_quadrature, chebyshev_gauss_integral, chebyshev_lobatto_integral

end
