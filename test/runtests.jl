using EllipseFit
using Test

@testset "EllipseFit" begin
    include("utils_tests.jl")
    include("ellipse_tests.jl")
    include("solver_tests.jl")
end