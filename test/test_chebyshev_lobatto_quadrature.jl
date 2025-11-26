@testset "cheb_lobatto_nodes" begin
    n = 8
    nodes = cheb_lobatto_nodes(n)
    @test length(nodes) == n + 1
    @test nodes[1] == 1.0
    @test nodes[end] == -1.0
    @test all(isapprox.(nodes, -reverse(nodes); atol=eps(Float64)))
    @test_throws ArgumentError cheb_lobatto_nodes(0)
end

@testset "cheb_lobatto_weights" begin
    n = 9
    weights = cheb_lobatto_weights(n)
    @test length(weights) == n + 1
    scale = 1 / n^2
    @test weights[1] == scale
    @test weights[end] == scale
    for k in 1:(n - 1)
        expected = isodd(k) ? 2 * scale : 0.0
        @test weights[k + 1] == expected
    end
    @test_throws ArgumentError cheb_lobatto_weights(0)
end

@testset "cheb_lobatto_quadrature" begin
    ns = (20, 40, 80)
    @testset "odd integrand vanishes" begin
        for n in ns
            @test isapprox(cheb_lobatto_quadrature(x -> x, n), 0.0; atol=eps(Float64))
        end
    end

    integrals = [
        (x -> 1.0, 2.0, "integral of 1"),
        (x -> x^2, 2/3, "integral of x^2"),
        (x -> sqrt(1 - x^2), pi / 2, "integral of sqrt"),
    ]
    bonus = (x -> exp(x), Base.MathConstants.e - inv(Base.MathConstants.e), "integral of exp")

    for (f, exact, label) in integrals
        @testset "$label" begin
            for n in ns
                approx = cheb_lobatto_quadrature(f, n)
                @test isapprox(approx, exact; rtol=0.995)
            end
        end
    end

    bonus_label = bonus[3]
    bonus_f = bonus[1]
    bonus_exact = bonus[2]
    @testset bonus_label begin
        for n in ns
            approx = cheb_lobatto_quadrature(bonus_f, n)
            @test isapprox(approx, bonus_exact; rtol=0.995)
        end
    end

    @test_throws ArgumentError cheb_lobatto_quadrature(x -> x, 0)
end
