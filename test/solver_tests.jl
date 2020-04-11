using EllipseFit
using LinearAlgebra
using Test


@testset "Test Gradient Desc" begin
    @testset "Quadratic Single Variable" begin
        dim = 1
        f = x -> (x .- 1).^2 .+ 2
        g = x -> 2 .* (x .- 1)
        a = 0.1
        xvals, fvals, gradnorms = graddesc(dim, f, g, alpha=a)

        @test isapprox(xvals[end], 1.0, atol=1e-6)
        @test isapprox(fvals[end], 2.0, atol=1e-6)
        @test isapprox(gradnorms[end], 0.0, atol=1e-6)
    end

    @testset "Quadratic Multi Variable" begin
        dim = 3
        f = x -> (2 * x[1] .- 1)^2 + 3 * (x[2] - 3)^2 + (x[3] + 1)^2 + 2
        g = x -> [2 * (2 * x[1] .- 1) * 2; 6 * (x[2] - 3); 2 * (x[3] + 1)]
        
        a = 0.01
        xvals, fvals, gradnorms = graddesc(dim, f, g, alpha=a)

        @test isapprox(xvals[:, end], [0.5, 3, -1], atol=1e-6)
    end
end


@testset "Test Gauss Newton" begin
    @testset "Quadratic" begin
        shp = (2,2)
        f = x -> [x[1]^2 - 1; x[2]^2 - 9]
        J = x -> [2*x[1] 0; 0 2*x[2]]
        a = 0.1
        xvals, fvals, gradnorms = gaussnewton(shp, f, J)
        @test isapprox(abs(xvals[1,end]), 1.0, atol=1e-6)
        @test isapprox(abs(xvals[2,end]), 3.0, atol=1e-6)
        @test isapprox(gradnorms[end], 0.0, atol=1e-6)
    end

    @testset "Range measurement problem" begin
        shp = (2,2)
        f = x -> [sqrt((x[1] - 2)^2 + (x[2] + 2)^2) - sqrt(10), sqrt((x[1] + 2)^2 + (x[2] - 3)^2) - sqrt(13)]
        J = x -> reshape([(x[1] - 2) / sqrt((x[1] - 2)^2 + (x[2] + 2)^2), (x[2] + 2) / sqrt((x[1] - 2)^2 + (x[2] + 2)^2),
                (x[1] + 2) / sqrt((x[1] + 2)^2 + (x[2] - 3)^2), (x[2] - 3) / sqrt((x[1] + 2)^2 + (x[2] - 3)^2)], (2,2))'

        xvals, fvals, gradnorms = gaussnewton(shp, f, J, xinit=[1, 0])
        @test isapprox(abs(xvals[1,end]), 1.0, atol=1e-6)
        @test isapprox(abs(xvals[2,end]), 1.0, atol=1e-6)
        @test isapprox(gradnorms[end], 0.0, atol=1e-6)
    end
end

@testset "Test Levenberg Marquardt" begin
    @testset "Quadratic" begin
        shp = (2,2)
        f = x -> [x[1]^2 - 1; x[2]^2 - 9]
        J = x -> [2*x[1] 0; 0 2*x[2]]
        a = 0.1
        xvals, fvals, gradnorms = levenbergmarquardt(shp, f, J)
        @test isapprox(abs(xvals[1,end]), 1.0, atol=1e-6)
        @test isapprox(abs(xvals[2,end]), 3.0, atol=1e-6)
        @test isapprox(gradnorms[end], 0.0, atol=1e-6)
    end

    @testset "Range measurement problem" begin
        shp = (2,2)
        f = x -> [sqrt((x[1] - 2)^2 + (x[2] + 2)^2) - sqrt(10), sqrt((x[1] + 2)^2 + (x[2] - 3)^2) - sqrt(13)]
        J = x -> reshape([(x[1] - 2) / sqrt((x[1] - 2)^2 + (x[2] + 2)^2), (x[2] + 2) / sqrt((x[1] - 2)^2 + (x[2] + 2)^2),
                (x[1] + 2) / sqrt((x[1] + 2)^2 + (x[2] - 3)^2), (x[2] - 3) / sqrt((x[1] + 2)^2 + (x[2] - 3)^2)], (2,2))'

        xvals, fvals, gradnorms = levenbergmarquardt(shp, f, J, xinit=[1, 0])
        @test isapprox(abs(xvals[1,end]), 1.0, atol=1e-6)
        @test isapprox(abs(xvals[2,end]), 1.0, atol=1e-6)
        @test isapprox(gradnorms[end], 0.0, atol=1e-6)
    end
end