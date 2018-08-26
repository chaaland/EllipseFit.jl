include("../utils/rmse.jl")


function gauss_newton(input_output_shape::Tuple{Int64,Int64}, f::Function, J::Function;
                      xinit=Inf, max_iters=1000, atol=1e-6)
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
        stop_criteria : the norm of the jacobian along the trajectory

    =#

    n = input_output_shape[1];
    m = input_output_shape[2];

    if any(isinf.(xinit))                 
        xinit = vec(randn(n));
    end
    
    xvals = hcat(xinit);
    xcurr = vec(xvals);
    fvals = vcat(f(xcurr));

    total_deriv = J(xcurr);
    stop_criteria = rmse(2*total_deriv' * f(xcurr));

    for i=1:max_iters
        if m == 1
            A = vcat(total_deriv);
            b = vcat(total_deriv .* xvals[:,i] .- fvals[:,i]);
            xcurr = A \ b;
        else
            A = vcat(total_deriv);
            b = vcat(total_deriv * xvals[:,i] - fvals[:,i]);
            xcurr = A \ b;
        end

        xvals = hcat(xvals, xcurr);
        fvals = hcat(fvals, f(xcurr));
        total_deriv = J(xcurr);
        stop_criteria = hcat(stop_criteria, rmse(2*total_deriv' * f(xcurr)));
        if rmse(2*total_deriv' * f(xcurr)) <= atol                   # From grad ||f(x)||^2
            break
        end
    end
    
    return xvals, fvals, vec(stop_criteria)
end