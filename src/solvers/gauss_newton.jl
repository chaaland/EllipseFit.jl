include("../utils/rmse.jl")

function gauss_newton(input_dim, f, J; xinit=Inf, max_iters=1000, atol=1e-6)
    #= Use the gauss-newton method to find extrema
    
    The Gauss-Newton method is used to approximately solve the non-linear least
    squares problem (NNLS). The method employs only first order information to 
    locate optima

    Args :
        input_dim : dimension of the input to the function to be minimized
        f : a function representing the residuals of the m nonlinear equations
        J : a function to evaluate the jacobian of 'f' 
        xinit : initial iterate for warm starting
        max_iters : the maximum number of newton steps to take
        atol : the absolute tolerance of the root mean square of the jacobian

    Returns :
        xvals : the trajectory of the gradient descent
        fvals : the value of the objective along the trajectory
        gradnorm : the norm of the jacobian along the trajectory

    =#

    if xinit == Inf
        xinit = randn((input_dim, 1))
    end

    xvals = xinit;
    fvals = f(xvals);
    total_deriv = J(xvals);
    gradnorm = norm(total_deriv)

    for i=1:max_iters
        gn_descent_direction = - (total_deriv' * total_deriv) \ (total_deriv * fvals)
        xcurr = xvals[:,i] + gn_descent_direction;
        xvals = hcat(xvals, xcurr);
        fvals = vcat(fvals, f(xcurr));
        total_deriv = J(xcurr);
        gradnorm = vcat(gradnorm, norm(g));
        if rmse(g) <= atol
            break
        end
    end
end