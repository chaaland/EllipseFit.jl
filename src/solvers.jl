
export Objective
export LeastSquares, OrthogonalDistance
export Solver
export NormalEquations, GradientDescent, GaussNewton, LevenbergMarquardt

# Objective subtypes
abstract type Objective end

struct LeastSquares <: Objective end
struct OrthogonalDistance <: Objective end

# Solver subtypes 
abstract type Solver end

struct NormalEquations <: Solver end

struct GradientDescent <: Solver
    alpha::Float64
    iterations::Int32 
    atol::Float64
    xinit

    function GradientDescent(alpha::Float64, iterations::Int32, atol::Float64, xinit=Inf)
        if alpha < 0
            error("Expected argument 'alpha' to be nonnegative")
        elseif iterations < 1
            error("Expected argument 'iterations' to be nonnegative")
        elseif atol <= 0
            error("Expected argument 'atol' to be positive")
        end
        return new(alpha, iterations, atol, xinit)
    end
end

struct GaussNewton <: Solver
    iterations::Int32 
    atol::Float64
    xinit

    function GaussNewton(iterations::Int32, atol::Float64, xinit=Inf)
        if iterations < 1
            error("Expected argument 'iterations' to be nonnegative")
        elseif atol <= 0
            error("Expected argument 'atol' to be positive")
        end
        return new(iterations, atol, xinit)
    end
end

struct LevenbergMarquardt <: Solver
    iterations::Int32
    atol::Float64
    xinit

    function LevenbergMarquardt(iterations::Int32, atol::Float64, xinit=Inf)
        if iterations < 1
            error("Expected argument 'iterations' to be nonnegative")
        elseif atol <= 0
            error("Expected argument 'atol' to be positive")
        end
        return new(iterations, atol, xinit)
    end
end

struct NewtonRaphson <: Solver
    iterations::Int32
    atol::Float64
    xinit

    function NewtonRaphson(iterations::Int32, atol::Float64, xinit=Inf)
        if iterations < 1
            error("Expected argument 'iterations' to be nonnegative")
        elseif atol <= 0
            error("Expected argument 'atol' to be positive")
        end
        return new(iterations, atol, xinit)
    end
end