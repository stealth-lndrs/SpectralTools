@testset "Poisson solver" begin
    f_fun(x, y) = -2 * pi^2 * sin(pi * x) * sin(pi * y)
    u_exact(x, y) = sin(pi * x) * sin(pi * y)
    for N in (10, 14, 18, 22)
        U, x, y = poisson_chebyshev_2d(f_fun, N)
        exact = [u_exact(xi, yi) for xi in x, yi in y]
        err = maximum(abs.(U - exact))
        @test err < 1e-4
    end
end
