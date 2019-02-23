@testset "Ellipse Creation" begin

    @test rmse(zeros(2)) == 0
    @test rmse(zeros(2,1)) ==0
    @test rmse(ones(4)) == 1
    @test rmse(ones(4), ones(4)) == 0
    @test rmse([1 -1], [-1 1]) == 2
    @test rmse([1; -1], [-1 1]) == 2

end
