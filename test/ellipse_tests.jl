using EllipseFit
using Test


@testset "Ellipse Quadratic Form Creation" begin

    S1 = [1 0; 0 1]
    S2 = [3 0; 0 2]
    c1 = vec([0 0])
    c2 = vec([1 -1])

    E1 = Ellipse(S1, c1)
    E2 = Ellipse(S2, c2)
    
    @test E1.quadform.S == S1
    @test E1.quadform.center == c1

    @test E1.conicform.A == 1
    @test E1.conicform.B == 0
    @test E1.conicform.C == 1
    @test E1.conicform.D == 0
    @test E1.conicform.E == 0

    @test E1.parametricform.center == c1
    @test E1.parametricform.semiaxis_lengths == vec([1 1])

    @test E2.quadform.S == S2
    @test E2.quadform.center == c2

    @test E2.conicform.A == 3
    @test E2.conicform.B == 0
    @test E2.conicform.C == 1
    # @test E2.conicform.D == 
    # @test E2.conicform.E == 
    # @test E2.conicform.F == 

    @test E2.parametricform.center == c2
    @test E2.parametricform.semiaxis_lengths == [3 2]
end
