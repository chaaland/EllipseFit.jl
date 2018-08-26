using PyPlot
using EllipseFit

# include("../src/utils/rotate_mat2d.jl");


function noisy_ellipse(semiaxis_lengths; center=[0 0], ccw_angle=0, numpoints=50)
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
    
    theta = 2 * pi * vec(rand(numpoints));
    onaxis_ellipse = vec(semiaxis_lengths) .* [cos.(theta) sin.(theta)]';
    epsilon = 0.3 * randn(2, numpoints);
    
    return  vec(center) .+ rotate_mat2d(ccw_angle) * onaxis_ellipse + epsilon;
end

X = noisy_ellipse([3 2], center = [1 -1], ccw_angle=-pi/3)

figure(figsize=(10,10));
scatter(X[1,:], X[2,:]);
xlim([-5,5]);
ylim([-5,5]);
grid(true);
savefig("../img/noisyellipse.png");