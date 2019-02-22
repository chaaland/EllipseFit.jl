using EllipseFit
using Test
using LinearAlgebra


@testset "vandermonde" begin

    @test size(vandermonde(vec([3 2 1]), 3)) == (3, 4)
    @test size(vandermonde(vec([3 2 1]), 0)) == (3, 1)
    @test_throws ErrorException vandermonde(vec([3 2 1]), -1)
    @test vandermonde(vec([2]), 2) == [1 2 4]
    @test vandermonde(vec([2 3]), 2) == [1 2 4; 1 3 9]
    @test vandermonde(vec([2 3]), 4) == [1 2 4 8 16; 1 3 9 27 81]
    @test vandermonde(vec([-2 -3]), 2) == [1 -2 4; 1 -3 9]

end