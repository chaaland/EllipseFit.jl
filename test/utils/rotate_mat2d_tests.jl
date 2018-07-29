using EllipsePlot

@testset "2d Rotations" begin

    @test size(EllipsePlot.rotate_mat2d(1)) == (2,2);
    @test EllipsePlot.rotate_mat2d(0) == eye(2,2);
    @test isapprox(EllipsePlot.rotate_mat2d(2 * pi) , eye(2,2));
    @test isapprox(EllipsePlot.rotate_mat2d(pi / 2), [0 -1; 1 0]);

    @test isapprox(EllipsePlot.rotate_mat2d(-pi / 2), 
                   EllipsePlot.rotate_mat2d(pi/2, ccw=false));

    @test isapprox(EllipsePlot.rotate_mat2d(pi/3) * EllipsePlot.rotate_mat2d(pi / 6),
                   EllipsePlot.rotate_mat2d(pi / 2));
                   
    @test isapprox(EllipsePlot.rotate_mat2d(1)' * EllipsePlot.rotate_mat2d(1), eye(2,2))
end

nothing