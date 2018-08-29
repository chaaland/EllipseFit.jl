include("../utils/utils.jl");


function levenbergmarquardt(input_output_shape::Tuple{Int64,Int64}, f::Function, J::Function; xinit=Inf, max_iters=1000, atol=1e-6)
    #= Implements the levenberg marquardt heuristic for finding roots of m nonlinear equations in n unknowns
    
    Args :
        input_output_shape : a tuple of the form (n, m) giving the number of variables
                             and the number of outputs respectively
        f : a function that takes 'input_dim' inputs and returns m outputs
        J : a function to evaluate the Jacobian of 'f'
        xinit : initial iterate for warm starting
        max_iters : the maximum number of iterations to perform 
        atol : the absolute tolerance of the root mean square of the Jacobian

    Returns :
        xvals : the trajectory of the gradient descent
        fvals : the value of the objective along the trajectory
        stop_criteria : root mean square of twice the transposed jacobian times the evaluation of the function
        lambdavals : the values of the penalty parameter for each iteration of the algo


        gradnorm : the norm of the gradient along the trajectory
        lambdavals : the values of the penalty parameter for each iteration of the algo
    =#

    LAMBDA_MAXIMUM = 1e10;
    STEPSZ_MINIMUM = 1e-5;

    n = input_output_shape[1];
    m = input_output_shape[2];

    if any(isinf.(xinit))                 
        xinit = vec(randn(n));
    end
    
    lambdavals = [1];
    xvals = hcat(xinit);
    xcurr = vec(xvals);
    fvals = vcat(f(xcurr));

    total_deriv = J(xcurr);
    stop_criteria = rmse(2 * total_deriv' * f(xcurr));

    for i in 1:max_iters
        while true
            if m == 1
                A = vcat(total_deriv, sqrt(lambdavals[i]) * Matrix{Float64}(I, n, n));
                b = vcat(total_deriv .* xvals[:,i] .- fvals[:,i], sqrt(lambdavals[i]) * xvals[:,i]);
                xcurr = A \ b;
            else
                A = vcat(total_deriv, sqrt(lambdavals[i]) * Matrix{Float64}(I, n, n));
                b = vcat(total_deriv * xvals[:,i] - fvals[:,i], sqrt(lambdavals[i]) * xvals[:,i]);
                xcurr = A \ b;
            end

            if norm(f(xcurr)) < norm(fvals[:,i])
                lambdavals = hcat(lambdavals, lambdavals[i] * 0.8);
                break
            elseif 2 * lambdavals[i] > LAMBDA_MAXIMUM
                lambdavals = hcat(lambdavals, lambdavals[i]);
                break;
            elseif rmse(xcurr - xvals[:,i]) < STEPSZ_MINIMUM
                lambdavals = hcat(lambdavals, lambdavals[i]);
                break;
            else
                lambdavals[i] = 2 * lambdavals[i];
            end
        end

        xvals = hcat(xvals, xcurr);
        fvals = hcat(fvals, f(xcurr));
        total_deriv = J(xcurr); 
        stop_criteria = hcat(stop_criteria, rmse(2 * total_deriv' * f(xcurr)));
        if stop_criteria[end] <= atol                   # From grad ||f(x)||^2
            break
        end
    end

    return xvals, fvals, vec(stop_criteria), vec(lambdavals)
end