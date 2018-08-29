using EllipsePlot

@testset "Ellipse coords from parametric form" begin
    A = eye(2,2);
    @test size(ellipse_from_parametric(A, numpoints=1000)) == (2, 1000);
    @test size(ellipse_from_parametric(A, numpoints=50)) == (2, 50);

end

nothing