include("rotate_mat2d.jl");


function parametric2ellipse_coords(semiaxis_lengths::Vector{T}; center=[0 0], ccw_angle=0, numpoints=1000) where T <: Real
    #= Helper for plotting a (possibly degenerate) 2D ellipse given semi-major/minor
    axes lengths, ellipse center coordinates and angle off the positive x axis
    
    Given the positve semidefinite matrix of a quadratic form specifying an ellipse
    in standard form x^T A x = 1, the center of the ellipse, the angle to rotate the
    ellipse wrt the positive x-axis, and the number of points desired for plotting, 
    an array containing the x, y coordinates of various points on the ellipse are returned
    
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

    center = vec(center);
    
    theta = vec(range(0, stop=2*pi, length=numpoints));
    onaxis_ellipse = semiaxis_lengths .* [cos.(theta) sin.(theta)]';
    
    return center .+ rotate_mat2d(ccw_angle) * onaxis_ellipse;
end