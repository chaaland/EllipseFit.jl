include("utils/elementwise_pseudoinvert.jl");

function quadform2ellipse_coords(psdMatrix; center=[0; 0], numpoints=1000)
    #= Helper for plotting a (possibly degenerate) 2D ellipse given its quadratic form
    
    Given the positve semidefinite matrix of a quadratic form specifying an ellipse
    in the form (x - c)^T A (x - c) = 1, the center of the ellipse 'c', and the number of points
    to plot, a matrix containing the x, y coordinates of various points on the ellipse
    are returned
    
    Args :
        psdMatrix : An array representing the quadratic form of an ellipse 
                    having the form 'x^T * psdMatrix * x = 1'
        center : The origin of the ellipse to be plotte. This is a 2x1 vector
        numpoints : The number of points to plot. Controls the granularity of the 
                    plotted figure
    
    Returns :
        A 2 x numpoints matrix where the x and y coordinates are in the 1st and 2nd row 
        respectively
    =#
    

    if psdMatrix != psdMatrix'
        error("Parameter 'psdMatrix' must be symmetric");
    elseif size(psdMatrix) != (2,2)
        error("Parameter 'psdMatrix' must be (2,2)")
    elseif size(vec(center))[0] != 2
        error("Parameter 'center' must be of size 2")
    end
    
    center = vec(center);
    if psdMatrix == zeros(2,2)
        return center * ones(1, numpoints)
    end

    f = eigfact(psdMatrix);
    V = f[:vectors];
    D = f[:values];
    
    D[D .== -0.0] = 0.0;
    negative_eigs = sum(D .< 0);
    if (negative_eigs > 0)
        error("Input matrix must be positive semidefinite");
    end
    
    theta = vec(linspace(0, 2*pi, numpoints));
    y = [cos.(theta) sin.(theta)]';         # trace out the unit circle
    invD = elementwise_pseudoinvert(D);
    X = V * (sqrt.(invD) .* y);             # stretch/rotate the circle into an ellipse using Sigma^{-1/2}
    
    return center .+ X
end
