using EllipsePlot

@testset "Ellipse coords from quadratic form" begin
    A = eye(2,2);
    @test size(ellipse_from_quadform(A, numpoints=1000)) == (2, 1000);
    @test size(ellipse_from_quadform(A, numpoints=50)) == (2, 50);

    X = ellipse_from_quadform(A, numpoints=1000); 
    @test diag(X' * A * X ) ≈ ones(1000);
    
    A = diagm([2, 4]);
    X = ellipse_from_quadform(A, numpoints=10); 
    @test diag(X' * A * X) ≈ ones(10);

    A = [4 1; 1 9];
    X = ellipse_from_quadform(A, numpoints=10); 
    @test diag(X' * A * X) ≈ ones(10);

    A = zeros(2,2);
    X = ellipse_from_quadform(A, center=[2 1], numpoints=10); 
    @test X ≈ zeros(2,10); 

    A = diagm([2, -1]);
    @test_throws ErrorException ellipse_from_quadform(A, numpoints=10); 

    c = [2 1];
    A = [4 1; 1 9];
    X = ellipse_from_quadform(A, center=c, numpoints=1000); 
    @test diag((X .- c')' * A * (X .- c')) ≈ ones(1000); 
 
end

nothing