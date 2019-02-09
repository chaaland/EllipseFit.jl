using EllipseFit
using Test
using LinearAlgebra


@testset "2d Rotations" begin

    @test size(rotation_mat(1)) == (2,2);
    @test rotation_mat(0) == Matrix{Float64}(I, 2, 2);
    @test isapprox(rotation_mat(2 * pi) , Matrix{Float64}(I, 2, 2));
    @test isapprox(rotation_mat(pi / 2), [0 -1; 1 0]);

    @test isapprox(rotation_mat(-pi / 2), 
                   rotation_mat(pi/2, ccw=false));

    @test isapprox(rotation_mat(pi/3) * rotation_mat(pi / 6),
                   rotation_mat(pi / 2));
                   
    @test isapprox(rotation_mat(1)' * rotation_mat(1), Matrix{Float64}(I, 2, 2))
end