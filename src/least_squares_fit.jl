
function least_squares_fit(X::Array{T,2}) where T <: Real
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

    N, n = size(X);

    if n != 2  
        error("Expected an array with second dimension 2")
    end

    x = vec(X[:,1]);
    y = vec(X[:,2]);
    one = ones(N,1);
    b = zeros(N, 1);

    A = hcat(x .^ 2, x .* y, y .^ 2, x, y);
    coeffs = A \ one;
    
    A = coeffs[1];
    B = coeffs[2];
    C = coeffs[3];
    D = coeffs[4];
    E = coeffs[5];
    
    return A, B, C, D, E, -1
end