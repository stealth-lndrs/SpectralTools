@testset "Chebyshev basics" begin
    x, D = cheb(16)
    @test length(x) == 17
    @test size(D) == (17, 17)

    f = x .^ 2
    df_numeric = D * f
    df_exact = 2 .* x
    @test maximum(abs.(df_numeric - df_exact)) < 1e-10

    weights = cheb_lobatto_weights(6)
    @test length(weights) == 7
    @test isapprox(weights[1], 1 / 36; atol = 0, rtol = 0)
    @test isapprox(weights[2], 2 / 36; atol = 0, rtol = 0)
    @test weights[3] == 0.0
end
