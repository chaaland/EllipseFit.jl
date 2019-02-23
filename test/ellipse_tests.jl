using EllipseFit
using Test


@testset "Ellipse Creation" begin

    S1 = [1 0; 0 1]
    c1 = vec([0 0])
    E1 = Ellipse(S1, c1)
    
    @test E1.quadform.S == S1
    @test E1.quadform.center == c1

    @test E1.conicform.A == 1
    @test E1.conicform.B == 0
    @test E1.conicform.C == 1
    @test E1.conicform.D == 0
    @test E1.conicform.E == 0

    @test E1.parametricform.center == c1
    @test E1.parametricform.semiaxis_lengths == vec([1 1])

end
