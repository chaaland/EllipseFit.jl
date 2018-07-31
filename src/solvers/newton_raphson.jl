include("../utils/rmse.jl")

function newton_raphson(input_dim, f, J, H; xinit=Inf, max_iters=10, atol=1e-6)
    #= Use the newton raphson method to find extrema
    
    The Newton-Raphson method takes into account second order information to 
    find locations where the function's gradient is zero (maxima, minima, 
    saddle points)

    Args :
        input_dim : dimension of the input to the function to be minimized
        f : a function that takes 'input_dim' inputs and returns a scalar
        J : a function to compute the jacobian of 'f' (i.e. the transpose of the gradient)
        H : hessian matrix of the function 'f' (i.e. the jacobian of the gradient)
        xinit : initial iterate for warm starting
        max_iters : the maximum number of newton steps to take
        atol : the absolute tolerance of the root mean square of the gradient

    Returns :
        xvals : the trajectory of the gradient descent
        fvals : the value of the objective along the trajectory
        gradnorm : the norm of the gradient along the trajectory

    =#

    if xinit == Inf
        xinit = randn((input_dim, 1))
    end

    xvals = xinit;
    fvals = f(xvals);
    g = J(xvals)';
    h = H(xvals);
    gradnorm = norm(g);

    for i=1:max_iters
        newton_step = - h \ g;
        xcurr = xvals[:,i] + newton_step;
        xvals = hcat(xvals, xcurr);
        fvals = vcat(fvals, f(xcurr));
        g = J(xcurr)';
        h = H(xcurr);
        gradnorm = vcat(gradnorm, norm(g));
        print(gradnorm)
        print(g)
        println(xvals)
        println(rmse(g))
        if rmse(g) <= atol
            break
        end
    end
    
    return xvals, fvals, gradnorm
end