using LinearAlgebra

include("utils.jl")
include("ellipse.jl")
include("solvers.jl")
include("solvers/levenbergmarquardt.jl")
include("solvers/gaussnewton.jl")

export EllipseModel, fit

mutable struct EllipseModel
    X::Array{Real,2}
    objective
    solver 

    function EllipseModel(
        X::Array{T,2},
        objective,
        solver,
    ) where {T <: Real}
        N, n = size(X)
        if n != 2
            error("Expected an array with second dimension 2")
        end
        new(X, objective, solver)
    end
end

function fit(model::EllipseModel)
    if model.objective == LeastSquares
       A, B, C, D, E = leastsquares(model.X, model.solver)
       return Ellipse(A, B, C, D, E)
    elseif model.objective == OrthogonalEuclideanDistance
        center, semiaxis_lengths, ccw_angle = orthogonaldist(model.X, model.solver)
        println(semiaxis_lengths)
        return Ellipse(semiaxis_lengths, center=center, ccw_angle=ccw_angle)
        # semiaxis_lengths::Array{T}; center=[0 0], ccw_angle=0
    else
        error("Unsupported objective type")
    end
end


function leastsquares(X::Array{T,2}, solver=NormalEquations) where T <: Real
    #= Fit an ellipse to data using least least squares fit

    Set up a least squares problem of the form ||A x - b||^2 where 
    A is the matrix of quadratic and linear terms, x is the vector of
    the coefficients of the ellipse represented in conic form and b is 
    a vector of all ones 

    Args :
        X : N x 2 matrix containing the data to be fit with an ellipse
            
    Returns :
        The coefficients of the ellipse in conic section form
    =#

    N, n = size(X)

    if n != 2  
        error("Expected an array with second dimension 2")
    end

    x = vec(X[:,1])
    y = vec(X[:,2])

    A = hcat(x .^ 2, x .* y, y .^ 2, x, y)
    b = ones(N)

    if typeof(solver) == NormalEquations
        coeffs = A \ b
        
        A = coeffs[1]
        B = coeffs[2]
        C = coeffs[3]
        D = coeffs[4]
        E = coeffs[5]
        
        return A, B, C, D, E
    elseif typeof(solver) == GradientDescent
        error("Unsupported solver GradientDescent")
    else
        error("Unsupported solver")
    end
end

function orthogonaldist(X::Array{T,2}, solver=LevenbergMarquardt) where T <: Real
    #= Fit an ellipse by minimizing the orthogonal distance of the points 

    Rather than least squares, the ellipse is fit so as to minimize the 
    perpendicular distance of all the measured points to an ellipse. This is
    an example of an errors in variables model which accounts for measurement
    error in both the independent and dependent variables

    Args :
        X : An N x 2 matrix containing the data to be fit with an ellipse

    Returns :
        center : a 2 vector containing x and y coords of ellipse center
        semiaxis_lengths : a 2 vector containing the semi axis lengths of the ellipse
        ccw_angle : scalar indicating the angle made wrt positive x axis

    =#

    N, n = size(X)
    if n != 2
        error("Expected array with second dimension 2")
    end

    n_params = 5 + N             # xcenter, ycenter, semi major, semi minor , angle + theta per data point
    n_equalities = 2 * N
    # println(solver.atol)
    if typeof(solver) == LevenbergMarquardt
        thetavals_lm, fvals_lm, gradnorm_lm, lambdavals_lm = levenbergmarquardt(
            (n_params, n_equalities), 
            z -> parametric_ellipse(z,X), 
            jacobian_ellipse; 
            xinit=solver.xinit,
            max_iters=solver.iterations,
            atol=solver.atol,
        )

        final_params = thetavals_lm[:,end]
        center = vec(final_params[1:2])
        semiaxis_lengths = vec(final_params[3:4])
        ccw_angle = final_params[5]

        return center, semiaxis_lengths, ccw_angle
    else
        error("Unsupported solver")
    end
end

function parametric_ellipse(x::Vector{T}, data) where T <: Real
    center = x[1:2]
    semiaxis_lengths = x[3:4]
    ccw_angle = x[5]
    theta = x[6:end]
    
    onaxis_ellipse = diagm(0 => vec(semiaxis_lengths)) * [cos.(theta)'; sin.(theta)']
    U = rotation_mat(ccw_angle)
    rotated_ellipse = U * onaxis_ellipse
    shifted_ellipse = center .+ rotated_ellipse

    return vec(shifted_ellipse) - vec(data')
end

function jacobian_ellipse(x::Vector{T}) where T <: Real
    center = vec(x[1:2])
    semiaxis_lengths = vec(x[3:4])
    ccw_angle = x[5]
    theta = vec(x[6:end])
    
    a = semiaxis_lengths[1]
    b = semiaxis_lengths[2]
    
    dxc = repeat([1, 0], length(theta), 1)
    dyc = repeat([0, 1], length(theta), 1)
    da = vec([cos.(ccw_angle)*cos.(theta) sin.(ccw_angle)*cos.(theta)]')
    db = vec([-sin.(ccw_angle)*sin.(theta) cos.(ccw_angle)*sin.(theta)]')
    dalpha = vec([-a*sin.(ccw_angle)*cos.(theta)-b*cos.(ccw_angle)*sin.(theta) a*cos.(ccw_angle)*cos.(theta)-b*sin.(ccw_angle)*sin.(theta)]')
    dtheta = zeros(2*length(theta), length(theta))
    
    for i = 1:length(theta) 
        theta_val = theta[i]
        val1 = -a * cos.(ccw_angle) * sin.(theta_val) - b * sin.(ccw_angle) * cos.(theta_val)
        val2 = -a * sin.(ccw_angle) * sin.(theta_val) + b * cos.(ccw_angle) * cos.(theta_val)
        dtheta[2 * i - 1:2 * i, i] = [val1 val2]
    end

    J = [dxc dyc da db dalpha dtheta]
end