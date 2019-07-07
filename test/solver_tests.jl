using EllipseFit
using LinearAlgebra
using Test


@testset "Test Gradient Desc" begin

    @testset "Unit circle" begin
        S = Matrix{Int32}(I, 2, 2)
        c = vec([0 0])
        E = Ellipse(S, c)

        @test E.quadform.S == S
        @test E.quadform.center == c

        @test E.conicform.A == 1
        @test E.conicform.B == 0
        @test E.conicform.C == 1
        @test E.conicform.D == 0
        @test E.conicform.E == 0

        @test E.parametricform.center == c
        @test E.parametricform.semiaxis_lengths == vec([1 1])
    end