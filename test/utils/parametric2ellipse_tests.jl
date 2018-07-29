using EllipsePlot

@testset "Ellipse coords from parametric form" begin
    A = eye(2,2);
    @test size(parametric2ellipse_coords(A, numpoints=1000)) == (2, 1000);
    @test size(parametric2ellipse_coords(A, numpoints=50)) == (2, 50);

end

nothing