include("../utils/rmse.jl")



#function levenberg_marquardt(input_output_shape::Tuple{Int64,Int64}, f::Function, J::Function; xinit=Inf, max_iters=1000, atol=1e-6)
function gauss_newton(input_output_shape::Tuple{Int64,Int64}, f::Function, J::Function; xinit=Inf, max_iters=1000, atol=1e-6)
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

    n = input_output_shape[1];
    m = input_output_shape[2];

    if any(isinf.(xinit))                 
        xinit = vec(randn(n));
    end
    
    xvals = hcat(xinit);
    xcurr = vec(xvals);
    fvals = vcat(f(xcurr));

    total_deriv = J(xcurr);
    gradnorm = norm(total_deriv);

    for i=1:max_iters
        if m == 1
            gn_step = fvals / total_deriv[0];
            xcurr = xvals[i] - gn_step;
        else
            A = total_deriv' * total_deriv; 
            b = total_deriv' * fvals[:,i];
            gn_step = A \ b;
            xcurr = xvals[:,i] - gn_step;
        end

        xvals = hcat(xvals, xcurr);
        fvals = vcat(fvals, f(xcurr));
        total_deriv = J(xcurr);
        gradnorm = vcat(gradnorm, norm(total_deriv));
        if rmse(total_deriv) <= atol
            break
        end
    end
end