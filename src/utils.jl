using Random;
using Statistics;

export rmse, rotation_mat, vandermonde, elementwise_pseudoinvert, noisyellipse

function rmse(x::AbstractArray, y::AbstractArray)
  #= Computes the root mean square error between x and y 

    Given two vectors x and y, return 
        sqrt(1/N * (sum((x_i - y_i)^2)))
    =#

    x = vec(x)
    y = vec(y)
    squared_diff = (x - y).^2
    mean_square_error = mean(squared_diff)

    return sqrt(mean_square_error)
end

function rmse(x::AbstractArray)
    #= Computes the root mean square error of x and y 

    Given two vectors x and y, return 
        sqrt(1/N * (sum((x_i - y_i)^2)))
    =#
    x = vec(x)
    squared_x = x.^2
    mean_square_error = mean(squared_x)

    return sqrt(mean_square_error)
end

function rotation_mat(angle::Real; ccw=true)
    #= Helper for rotating points in 2d space
    
    Builds the rotation matrix specified by the angle and the direction of rotation
    
    Args :
        angle : angle in radians of rotation
        ccw : whether the angle is measured counter clockwise wrt the positive x axis
    
    Returns :
        A 2x2 rotation matrix
    =#

    rotate_mat = zeros(2,2);
    if ccw
        rotate_mat = [cos(angle) -sin(angle);
                      sin(angle) cos(angle)];
    else
        rotate_mat = [cos(angle) sin(angle);
                      -sin(angle) cos(angle)];
    end
    
    return rotate_mat;
end

function vandermonde(x::Vector{T}, degree::Int) where T <: Real
    #= Generates a vandermonde matrix using the entries of x

    Given a vector x, the successive powers of x up to and including n
    are computed and stored as rows of a matrix

    Args :
        x : a vector of values
        n : the degree of the polynomial
    Returns :
        A matrix of dimension length(x) by n+1 where a row is given by
            [1 a a^2 ... a^{n-1} a^n]
    =#

    if degree < 0
        error("Polynomial degree is expected to be nonnegative")
    end 
    N = length(x)
    A = zeros(N, degree + 1)
    col = ones(N)

    for i = 0:degree
        A[:, i+1] = col;
        col = col .* x;
    end

    return A
end

function elementwise_pseudoinvert(v::AbstractArray, tol=1e-10) 
    #= Helper for elementwise inversion of a matrix v
    
    Inverts the non-zero elements of v and keeps the 0 elements. 
    Elements within tol of 0 are treated as zero
    
    Args :
        v : an array to take the reciprocal of
        tol : tolerance for setting an element to zero
    
    Returns :
        Array with the elements inverted except the zeros
    =#
    
    m = maximum(abs.(v));
    if m == 0
        return v
    end
    v = v ./ m;
    reciprocal = 1 ./ v;
    reciprocal[abs.(reciprocal) .>= 1/tol] .= 0;   

    return reciprocal / m;
end

function noisyellipse(semiaxis_lengths; center=[0 0], ccw_angle=0, numpoints=50)
    #= Example function for plotting data roughly corresponding to an ellipse plus noise

    Helper function used to generate noisy ellipse data with sample points drawn from 
    a uniform radial distribution

    Args :
        semiaxis_lengths : array of the form [a b] where a is half the length of the axis aligned
                      ellipse along the x-axis and b is half the length along the y-axis (before
                      rotation)
        center : array of the x and y coordinates of the center of the ellipse
        ccw_angle : The counter clockwise angle (in rad) to rotate the ellipse wrt
                    the positive x-axis
        num_points : an integer indicating the number of xy pairs on the ellipse to return    
    Returns :
        A 2 x numpoints array where the x and y coordinates are in the 1st and 2nd row 
        respectively
    =# 

    if size(vec(center))[1] != 2
        error("Parameter 'center' must be of size 2");
    end
    
    if size(vec(semiaxis_lengths))[1] != 2
        error("Parameter 'semiaxis_lengths' must be of size 2");
    end
    
    theta = 2 * pi * vec(rand(numpoints))
    onaxis_ellipse = vec(semiaxis_lengths) .* [cos.(theta) sin.(theta)]'
    epsilon = 0.3 * randn(2, numpoints)
    noisy_ellipse_coords = vec(center)' .+ (rotation_mat(ccw_angle) * onaxis_ellipse)' + epsilon'

    return noisy_ellipse_coords
end