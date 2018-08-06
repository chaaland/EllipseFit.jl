include("../utils/rmse.jl");


function levenberg_marquardt(input_output_shape::Tuple{Int64,Int64}, f::Function, J::Function; xinit=Inf, max_iters=1000, atol=1e-6)
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
        gradnorm : the norm of the gradient along the trajectory
        lambdavals : the values of the penalty parameter for each iteration of the algo

    =#

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
    gradnorm = rmse(2*total_deriv' * f(xcurr));

    for i in 1:max_iters
        while true
            if m == 1
                A = total_deriv[1]^2 + lambdavals[i];          
                b = total_deriv[1] * fvals[i];
                lm_step = b / A ;
                xcurr = xvals[i] - lm_step;
            else
                A = total_deriv' * total_deriv + lambdavals[i] * eye(n, n); 
                b = total_deriv' * fvals[:,i];
                lm_step = A \ b;
                xcurr = xvals[:,i] - lm_step;
            end

            if sum(f(xcurr).^2) > sum(fvals[:,i].^2)
                lambdavals[i] = 2 * lambdavals[i];
            else
                lambdavals = hcat(lambdavals, lambdavals[i] * 0.8);
                break
            end
        end

        xvals = hcat(xvals, xcurr);
        fvals = hcat(fvals, f(xcurr));
        total_deriv = J(xcurr); 
        gradnorm = hcat(gradnorm, rmse(2*total_deriv' * f(xcurr)));
        if rmse(2*total_deriv' * f(xcurr)) <= atol
            break
        end
    end
    
    return xvals, fvals, gradnorm, lambdavals
end