function rotate_mat2d(angle; ccw=true)
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
