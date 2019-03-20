using EllipseFit
using LinearAlgebra
using Test


@testset "Ellipse Quadratic Form Creation" begin

    @testset "Unit circle" begin
        S = Matrix{Int32}(I, 2, 2)
        c = vec([0 0])
        E = Ellipse(S, c)

        @test E.quadform.S == S
        @test E.quadform.center == c

        @test E.conicform.A == 1
        @test E.conicform.B == 0
        @test E.conicform.C == 1
        @test E.conicform.D == 0
        @test E.conicform.E == 0

        @test E.parametricform.center == c
        @test E.parametricform.semiaxis_lengths == vec([1 1])
    end

    @testset "Horizontal Ellipse Origin" begin
        a = 3
        b = 2
        S = diagm(0 => vec([1/a^2, 1/b^2]))
        c = vec([0 0])
        E = Ellipse(S, c)

        @test E.quadform.S == S
        @test E.quadform.center == c

        @test E.conicform.A == 1/a^2
        @test E.conicform.B == 0
        @test E.conicform.C == 1/b^2
        @test E.conicform.D == 0
        @test E.conicform.E == 0

        @test E.parametricform.center == c
        @test isapprox(E.parametricform.semiaxis_lengths, vec([a, b]))
        @test isapprox(E.parametricform.ccw_angle % pi, 0)
    end

    @testset "Vertical Ellipse Origin" begin
        a = 2
        b = 3
        S = diagm(0 => vec([1/a^2, 1/b^2]))
        c = vec([0 0])
        E = Ellipse(S, c)

        @test E.quadform.S == S
        @test E.quadform.center == c

        @test E.conicform.A == 1/a^2
        @test E.conicform.B == 0
        @test E.conicform.C == 1/b^2
        @test E.conicform.D == 0
        @test E.conicform.E == 0

        @test E.parametricform.center == c
        @test isapprox(E.parametricform.semiaxis_lengths, vec([b, a]))
        @test isapprox(E.parametricform.ccw_angle % pi, pi / 2)
    end

    @testset "Horizontal Ellipse Off Center" begin
        a = 3
        b = 2
        S = diagm(0 => vec([1/a^2, 1/b^2]))
        c = vec([1 -1])
        E = Ellipse(S, c)

        @test E.quadform.S == S
        @test E.quadform.center == c

        negF = 1 - 1/a^2 - 1/b^2
        @test E.conicform.A == (1 / a^2) / negF
        @test E.conicform.B == 0
        @test E.conicform.C == (1 / b^2) / negF
        @test E.conicform.D == (-2 * 1/a^2) / negF
        @test E.conicform.E == (-2 * -1/b^2) / negF

        @test E.parametricform.center == c
        @test isapprox(E.parametricform.semiaxis_lengths, [a, b])
        @test isapprox(E.parametricform.ccw_angle % pi, 0)
    end

    @testset "Vertical Ellipse Off Center" begin
        a = 2
        b = 3
        S = diagm(0 => vec([1/a^2, 1/b^2]))
        c = vec([1 -1])
        E = Ellipse(S, c)

        @test E.quadform.S == S
        @test E.quadform.center == c

        negF = 1 - 1/a^2 - 1/b^2
        @test E.conicform.A == (1 / a^2) / negF
        @test E.conicform.B == 0
        @test E.conicform.C == (1 / b^2) / negF
        @test E.conicform.D == (-2 * 1/a^2) / negF
        @test E.conicform.E == (-2 * -1/b^2) / negF

        @test E.parametricform.center == c
        @test isapprox(E.parametricform.semiaxis_lengths, [b, a])
        @test isapprox(E.parametricform.ccw_angle % pi, pi/2)
    end

    @testset "CCW Rotation 45 deg" begin
        a = 3
        b = 2
        ccw_angle = pi / 4
        rot_mat = rotation_mat(ccw_angle)
        S = rot_mat * diagm(0 => [1/a^2; 1/b^2]) * rot_mat'
        c = vec([0 0])

        E = Ellipse(S, c)
        @test E.quadform.S == S
        @test E.quadform.center == c

        @test E.conicform.A == S[1,1]
        @test E.conicform.B == (S[2,1] + S[1,2])
        @test E.conicform.C == S[2,2]
        @test E.conicform.D == 0
        @test E.conicform.E == 0

        @test E.parametricform.center == c
        @test isapprox(E.parametricform.semiaxis_lengths, vec([a b]))
        @test isapprox((E.parametricform.ccw_angle + 2*pi)% pi, ccw_angle)
    end

    @testset "CCW Rotation -45 deg" begin
        a = 3
        b = 2
        ccw_angle = -pi / 4
        rot_mat = rotation_mat(ccw_angle)
        S = rot_mat * diagm(0 => [1/a^2; 1/b^2]) * rot_mat'
        c = vec([0 0])

        E = Ellipse(S, c)
        @test E.quadform.S == S
        @test E.quadform.center == c

        @test E.conicform.A == S[1,1]
        @test E.conicform.B == (S[2,1] + S[1,2])
        @test E.conicform.C == S[2,2]
        @test E.conicform.D == 0
        @test E.conicform.E == 0

        @test E.parametricform.center == c
        @test isapprox(E.parametricform.semiaxis_lengths, vec([a b]))
        @test isapprox((E.parametricform.ccw_angle + 2*pi) % pi, (ccw_angle + 2*pi) % pi)
    end
end

@testset "Ellipse Parametric Form Creation" begin
    @testset "Unit Circle" begin 
        semiaxes = vec([1 1])
        E = Ellipse(semiaxes, center=[0 0]) 

        @test E.conicform.A == 1
        @test E.conicform.B == 0
        @test E.conicform.C == 1
        @test E.conicform.D == 0
        @test E.conicform.E == 0

        @test E.quadform.center == vec([0 0])
        @test E.quadform.S == [1 0; 0 1]
    end

    @testset "Unit Circle Off Center" begin 
        semiaxes = vec([1 1])
        c = vec([1 2])
        E = Ellipse(semiaxes, center=c) 

        negF = 1 - 1 - 2^2
        @test E.conicform.A == 1 / negF
        @test E.conicform.B == 0
        @test E.conicform.C == 1 / negF
        @test E.conicform.D == -2 * c[1] / negF
        @test E.conicform.E == -2 * c[2] / negF

        @test E.quadform.center == c
        @test E.quadform.S == Matrix(I,2,2)
    end

    @testset "Horizontal Ellipse" begin 
        a = 3
        b = 2
        semiaxes = vec([a b])
        c = vec([0 0])
        E = Ellipse(semiaxes, center=c) 

        @test E.conicform.A == 1/a^2
        @test E.conicform.B == 0
        @test E.conicform.C == 1/b^2
        @test E.conicform.D == 0
        @test E.conicform.E == 0

        @test E.quadform.center == c
        @test E.quadform.S == diagm(0 => [1/a^2; 1/b^2])
    end

    @testset "Vertical Ellipse 1" begin 
        a = 2
        b = 3
        semiaxes = vec([a b])
        c = vec([0 0])
        E = Ellipse(semiaxes, center=c) 

        @test E.conicform.A == 1/a^2
        @test E.conicform.B == 0
        @test E.conicform.C == 1/b^2
        @test E.conicform.D == 0
        @test E.conicform.E == 0

        @test E.quadform.center == c
        @test isapprox(E.quadform.S, diagm(0 => [1/a^2; 1/b^2]))
    end

    @testset "Vertical Ellipse 2" begin 
        a = 3
        b = 2
        semiaxes = vec([a b])
        c = vec([0 0])
        angle = pi / 2
        E = Ellipse(semiaxes, center=c, ccw_angle=angle)

        @test E.conicform.A == 1/b^2
        @test isapprox(E.conicform.B, 0., atol = 1e-5)
        @test E.conicform.C == 1/a^2
        @test E.conicform.D == 0
        @test E.conicform.E == 0

        @test E.quadform.center == c
        @test isapprox(E.quadform.S, diagm(0 => [1/b^2; 1/a^2]))
    end

    @testset "CCW Rotation 45 deg" begin 
        a = 2
        b = 3
        semiaxes = vec([a b])
        c = vec([0 0])
        angle = pi / 4
        E = Ellipse(semiaxes, center=c, ccw_angle=angle)

        S = rotation_mat(angle) * diagm(0 => 1 ./semiaxes.^2) * rotation_mat(angle)'

        @test E.conicform.A == S[1,1]
        @test isapprox(E.conicform.B, (S[1,2] + S[2,1]), atol = 1e-5)
        @test E.conicform.C == S[2,2]
        @test E.conicform.D == 0
        @test E.conicform.E == 0

        @test E.quadform.center == c
        @test isapprox(E.quadform.S, S)
    end

    @testset "CCW Rotation -45 deg" begin 
        a = 2
        b = 3
        semiaxes = vec([a b])
        c = vec([0 0])
        angle = -pi / 4
        E = Ellipse(semiaxes, center=c, ccw_angle=angle)

        S = rotation_mat(angle) * diagm(0 => 1 ./semiaxes.^2) * rotation_mat(angle)'

        @test E.conicform.A == S[1,1]
        @test isapprox(E.conicform.B, (S[1,2] + S[2,1]), atol = 1e-5)
        @test E.conicform.C == S[2,2]
        @test E.conicform.D == 0
        @test E.conicform.E == 0

        @test E.quadform.center == c
        @test isapprox(E.quadform.S, S)
    end

end

@testset "Ellipse Concic Form Creation" begin
    @testset "Unit Circle" begin 
        A = 1
        B = 0
        C = 1
        D = 0
        E = 0

        E = Ellipse(A, B, C, D, E)

        @test E.parametricform.center == vec([0 0])
        @test isapprox(E.parametricform.semiaxis_lengths, vec([1 1]))

        @test E.quadform.center == vec([0 0])
        @test E.quadform.S == [1 0; 0 1]
    end

    @testset "Horizontal Ellipse" begin 
        A = 1/4
        B = 0
        C = 1
        D = 0
        E = 0

        ellipse = Ellipse(A, B, C, D, E)

        @test ellipse.parametricform.center == vec([0 0])
        @test isapprox(ellipse.parametricform.semiaxis_lengths, vec([2 1]))
        @test isapprox((ellipse.parametricform.ccw_angle + 2*pi) % pi, 0)

        @test ellipse.quadform.center == vec([0 0])
        @test ellipse.quadform.S == [A B/2; B/2 C]
    end

    @testset "Vertical Ellipse" begin 
        A = 1
        B = 0
        C = 1/4
        D = 0
        E = 0

        ellipse = Ellipse(A, B, C, D, E)

        @test ellipse.parametricform.center == vec([0 0])
        @test isapprox(ellipse.parametricform.semiaxis_lengths, vec([2 1]))
        @test isapprox((ellipse.parametricform.ccw_angle + 2*pi) % pi, pi/2)

        @test ellipse.quadform.center == vec([0 0])
        @test ellipse.quadform.S == [A B/2; B/2 C]
    end

    @testset "Horizontal off center" begin 
        A = 1/4 / (-1/4)
        B = 0
        C = 1 / (-1/4)
        D = -2 * 1/4 / (-1/4)
        E = -2 * -1 / (-1/4)

        ellipse = Ellipse(A, B, C, D, E)

        @test ellipse.parametricform.center == vec([1 -1])
        @test isapprox(ellipse.parametricform.semiaxis_lengths, vec([2 1]))
        @test isapprox((ellipse.parametricform.ccw_angle + 2*pi) % pi, 0)

        @test ellipse.quadform.center == vec([1 -1])
        @test ellipse.quadform.S == diagm(0 => vec([1/4 1]))
    end

    @testset "CCW Rotation 45 deg" begin 
        angle = pi / 4
        a = 3
        b = 2
        semiaxes = vec([a b])
        S = rotation_mat(angle) * diagm(0 => 1 ./semiaxes.^2) * rotation_mat(angle)'
        A = S[1,1]
        B = (S[1,2] + S[2,1])
        C = S[2,2]
        D = 0
        E = 0

        ellipse = Ellipse(A, B, C, D, E)

        @test ellipse.parametricform.center == vec([0 0])
        @test isapprox(ellipse.parametricform.semiaxis_lengths, vec([a b]))
        @test isapprox((ellipse.parametricform.ccw_angle + 2*pi) % pi, angle)

        @test ellipse.quadform.center == vec([0 0])
        @test ellipse.quadform.S == [A B/2; B/2 C]
    end

    @testset "CCW Rotation -45 deg" begin 
        angle = -pi / 4
        a = 3
        b = 2
        semiaxes = vec([a b])
        S = rotation_mat(angle) * diagm(0 => 1 ./semiaxes.^2) * rotation_mat(angle)'
        A = S[1,1]
        B = (S[1,2] + S[2,1])
        C = S[2,2]
        D = 0
        E = 0

        ellipse = Ellipse(A, B, C, D, E)

        @test ellipse.parametricform.center == vec([0 0])
        @test isapprox(ellipse.parametricform.semiaxis_lengths, vec([a b]))
        @test isapprox((ellipse.parametricform.ccw_angle + 2*pi) % pi, (angle + 2*pi) % pi)

        @test ellipse.quadform.center == vec([0 0])
        @test ellipse.quadform.S == [A B/2; B/2 C]
    end
end