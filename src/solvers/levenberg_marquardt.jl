include("../utils/rmse.jl");


function levenberg_marquardt(input_output_shape, f, J; xinit=Inf, max_iters=1000, atol=1e-6)
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

    if isinf(xinit)
        xinit = vec(randn(n));
    end
    
    lambdavals = hcat(1);
    xvals = hcat(xinit);
    xcurr = xinit;
    fvals = hcat(f(xvals));
    total_deriv = J(xvals);
    gradnorm = norm(total_deriv);

    for i in 1:max_iters
        while true
            if m == 1
                A = total_deriv^2 + lambdavals[i];          
                b = total_deriv * fvals[i];
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
        gradnorm = vcat(gradnorm, norm(total_deriv));
        if rmse(total_deriv) <= atol
            break
        end
    end
    
    return xvals, fvals, gradnorm, lambdavals
end