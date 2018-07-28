include("../utils/rmse.jl")

function grad_desc(input_dim, f, grad; alpha=0.1, xinit=Inf, max_iters=1000, atol=1e-6)
    #= Perform gradient descent to find minimum of a function
    
    Args :
        input_dim : dimension of the input to the function to be minimized
        f : a function that takes 'input_dim' inputs and returns a scalar
        grad : gradient of the function 'f'
        alpha : the learning rate determining the step size
        max_iters : the maximum number of gradient steps to take
        atol : the absolute tolerance of the root mean square of the gradient

    Returns :
        xvals : the trajectory of the gradient descent
        fvals : the value of the objective along the trajectory
        gradnorm : the norm of the gradient along the trajectory

    =#

    if xinit == Inf
        xvals = randn((input_dim, 1));
    end

    xvals = xinit;
    fvals = f(xvals);
    g = grad(xvals);
    gradnorm = norm(g);

    for i in 1:max_iters
        xcurr = xvals[:,i] - alpha * g;
        xvals = hcat(xvals, xcurr);
        fvals = vcat(fvals, f(xcurr));
        g = grad(xcurr); 
        gradnorm = vcat(gradnorm, norm(g));
        if rmse(g) <= atol
            break
        end
    end

    return xvals, fvals, gradnorm
end