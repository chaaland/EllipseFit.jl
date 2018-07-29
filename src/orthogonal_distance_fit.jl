
function orthogonal_dist_fit(X)
    #= Fit an ellipse by minimizing the orthogonal distance of the points 

    Rather than least squares, the ellipse is fit so as to minimize the 
    perpendicular distance of all the measured points an ellipse. This is
    an example of an errors in variables model which accounts for measurement
    error in both the independent and dependent variables

    Args :
        X : 2 x N or N x 2 matrix containing the data to be fit with an ellipse

    Returns :

    =#
    
end