export Objective
export LeastSquares, OrthogonalEuclideanDistance
export Solver
export NormalEquations, GradientDescent, GaussNewton, LevenbergMarquardt

# Objective subtypes
abstract type Objective end

struct LeastSquares <: Objective end
struct OrthogonalEuclideanDistance <: Objective end

# Solver subtypes 
abstract type Solver end

struct NormalEquations <: Solver end

struct GradientDescent <: Solver
    alpha::Float64
    iterations::Int32 
    atol::Float64
    xinit::Union{Array{Real}, Real}

    function GradientDescent(alpha::Float64, iterations::Int32, atol::Float64, xinit=Inf::Union{Array{Real}, Real})
        if alpha < 0
            error("Expected argument 'alpha' to be nonnegative")
        elseif iterations < 1
            error("Expected argument 'iterations' to be positive")
        elseif atol <= 0
            error("Expected argument 'atol' to be positive")
        end
        return new(alpha, iterations, atol, xinit)
    end
end

struct GaussNewton <: Solver
    iterations::Int32 
    atol::Float64
    xinit::Union{Array{Real}, Real}

    function GaussNewton(iterations::Int32, atol::Float64, xinit=Inf::Union{Array{Real}, Real})
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
    xinit::Union{Array{Real}, Real}

    function LevenbergMarquardt(iterations::Int32, atol::Float64, xinit=Inf::Union{Array{Real}, Real})
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
    xinit::Union{Array{Real}, Real}

    function NewtonRaphson(iterations::Int32, atol::Float64, xinit=Inf::Union{Array{Real}, Real})
        if iterations < 1
            error("Expected argument 'iterations' to be nonnegative")
        elseif atol <= 0
            error("Expected argument 'atol' to be positive")
        end
        return new(iterations, atol, xinit)
    end
end