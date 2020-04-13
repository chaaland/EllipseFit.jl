using EllipseFit
using LinearAlgebra


@testset "root mean square test" begin

    @test rmse(zeros(2)) == 0
    @test rmse(zeros(2,1)) ==0
    @test rmse(ones(4)) == 1
    @test rmse(ones(4), ones(4)) == 0
    @test rmse([1 -1], [-1 1]) == 2
    @test rmse([1; -1], [-1 1]) == 2

end

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

@testset "vandermonde" begin

    @test size(vandermonde(vec([3 2 1]), 3)) == (3, 4)
    @test size(vandermonde(vec([3 2 1]), 0)) == (3, 1)
    @test_throws ErrorException vandermonde(vec([3 2 1]), -1)
    @test vandermonde(vec([2]), 2) == [1 2 4]
    @test vandermonde(vec([2 3]), 2) == [1 2 4; 1 3 9]
    @test vandermonde(vec([2 3]), 4) == [1 2 4 8 16; 1 3 9 27 81]
    @test vandermonde(vec([-2 -3]), 2) == [1 -2 4; 1 -3 9]

end

@testset "Elementwise Pseudo-inversion" begin

    A = [1 0; -1 2.5]
    B = [1 0; -1 2.5; 0 0]
    @test size(elementwise_pseudoinvert(A)) == size(A)
    @test size(elementwise_pseudoinvert(B)) == size(B)
    @test elementwise_pseudoinvert(zeros(2,2)) == zeros(2,2)
    @test isapprox(elementwise_pseudoinvert(A), [1 0; -1 0.4])
    @test isapprox(elementwise_pseudoinvert(B), [1 0; -1 0.4; 0 0])

end