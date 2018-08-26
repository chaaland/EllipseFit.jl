using LinearAlgebra; 

include("elementwise_pseudoinvert.jl")
include("rotate_mat2d.jl")


function quad2conic(S::AbstractArray; center=[0 0])
    #= Helper for converting quadratic form ellipse into a standard conic form
    
    Given an ellipse as a quadratic form

            (x - c)^T S (x - c) = 1

    convert it to the conic section form 

        A * x^2 + B * x * y + C * y^2 + D * x + E * y + F = 0

    This is done by completing the square leading to the following equivalence
            x^T * Q * x + b^T * x + c = 0
    
    Args :
        S : PSD matrix representing the ellipse in the form (x - c)^T S (x - c)
        center : The origin of the ellipse to be plotte. This is a 2x1 vector

    Returns :
        A : Coefficient of x^2
        B : Coefficient of x * y 
        C : Coefficient of y^2
        D : Coefficient of x
        E : Coefficient of y 
        F : Bias
    =# 

    A = S[1,1];
    B = S[1,2] + S[2,1];
    C = S[2,2];

    b = -2 * center' * S;
    D = b[1];
    E = b[2];
    F = center' * S * center - 1;

    return A, B, C, D, E, F
end

function conic2quad(A, B, C, D, E, F)
    #= Helper for converting conic form ellipse into a standard quadratic form
    
    Given an ellipse as a conic section of the form

        A * x^2 + B * x * y + C * y^2 + D * x + E * y + F = 0

    convert it to the form 

        (x - c)^T S (x - c) = 1

    This is done by completing the square leading to the following equivalence
            x^T * Q * x + b^T * x + c = 0
        =>  (x + 0.5 * inv(Q) * b)^T Q (x + 0.5 * inv(Q) * b)) + c - 0.25 * b^T * inv(A) * b = 0
    
    Args :
        A : Coefficient of x^2
        B : Coefficient of x * y 
        C : Coefficient of y^2
        D : Coefficient of x
        E : Coefficient of y 
        F : Bias

    Returns :
        S : PSD matrix representing the ellipse in the form (x - c)^T S (x - c)
        center : The origin of the ellipse to be plotte. This is a 2x1 vector

    =#

    if B^2 - 4*A*C >= 0
        error("Discriminant is non-negative. Input does not denote an ellipse");
    end

    Q = [A B/2; B/2 C];
    b = [D; E];

    beta = Q \ b;
    c = 0.25 * b' * beta - F;

    S = Q / c;
    center = -0.5 * beta;

    return S, center
end

function parametric2quad(semiaxis_lengths; center=[0 0], ccw_angle=0)
    #= Helper for converting from paramteric to quadratic form of ellipse

    Given an ellipse in parametric form

        [x_c y_c] + rot_mat2d(ccw_angle) * [a*cos(theta) b*sin(theta)]

    covnert it into a quadratic of the form

        (x - c)^T S (x - c) = 1

    Args :
        semiaxis_lengths : array of the form [a b] where a is half the length
                           of the axis aligned ellipse along the x-axis and b
                           is half the length along the y-axis (before rotation)
        center : array of the x and y coordinates of the center of the ellipse
        ccw_angle : The counter clockwise angle (in rad) to rotate the ellipse wrt
                    the positive x-axis

    Returns :
       S : PSD matrix representing the ellipse in the form (x - c)^T S (x - c)
       center : The origin of the ellipse to be plotte. This is a 2x1 vector 
    =#

    if size(vec(center))[1] != 2
        error("Parameter 'center' must be of size 2");
    end

    sqrtD = diagm(0 => vec(semiaxis_lengths));
    invD = elementwise_pseudoinvert(sqrtD .^ 2);
    V = rotate_mat2d(ccw_angle);
    S = V * invD * V';

    return S, vec(center)
end

function quad2parametric(S; center=[0 0])
    #= Helper for converting from quadratic form to parameteric form of ellipse

    Given an ellipse as a quadratic form 

            (x - c)^T S (x - c) = 1

    convert it to the parametric form

            [x_c y_c] + rot_mat2d(ccw_angle) * [a*cos(theta) b*sin(theta)]

    Args :
        S : PSD matrix representing the ellipse in the form (x - c)^T S (x - c)
        center : The origin of the ellipse to be plotte. This is a 2x1 vector 

    Returns :
        semiaxis_lengths : array of the form [a b] where a is half the length
                           of the axis aligned ellipse along the x-axis and b
                           is half the length along the y-axis (before rotation)
        center : array of the x and y coordinates of the center of the ellipse
        ccw_angle : The counter clockwise angle (in rad) to rotate the ellipse wrt
                    the positive x-axis
    =#

    f = eigen(S);
    V = f.vectors;
    D = f.values;

    semiaxis_lengths = sqrt(elementwise_pseudoinvert(D));
    ccw_angle = acos(V[0,0]);

    return vec(semiaxis_lengths), vec(center), ccw_angle
end

function parametric2conic(semiaxis_lengths; center=[0 0], ccw_angle=0)
    #= Helper for converting from parametric form to standard form of ellipse

    Given an ellipse in parametric form

            [x_c y_c] + rot_mat2d(ccw_angle) * [a*cos(theta) b*sin(theta)]

    convert it to the standard conic section form

            A * x^2 + B * x * y + C * y^2 + D * x + E * y + F = 0

    Args :
        semiaxis_lengths : array of the form [a b] where a is half the length
                           of the axis aligned ellipse along the x-axis and b
                           is half the length along the y-axis (before rotation)
        center : array of the x and y coordinates of the center of the ellipse
        ccw_angle : The counter clockwise angle (in rad) to rotate the ellipse wrt
                    the positive x-axis

    Returns :
        A : Coefficient of x^2
        B : Coefficient of x * y 
        C : Coefficient of y^2
        D : Coefficient of x
        E : Coefficient of y 
        F : Bias        
    =#

    S, center = parametric2quad(semiaxis_lengths, center, ccw_angle);
    return quad2conic(S, center);
end

function conic2parametric(A, B, C, D, E, F)
    #= Helper for converting from conic form to parametric form of ellipse

    Given an ellipse in standard conic section form

            A * x^2 + B * x * y + C * y^2 + D * x + E * y + F = 0

    convert it to parametric form

            [x_c y_c] + rot_mat2d(ccw_angle) * [a*cos(theta) b*sin(theta)]

    Args :
        A : Coefficient of x^2
        B : Coefficient of x * y 
        C : Coefficient of y^2
        D : Coefficient of x
        E : Coefficient of y 
        F : Bias 

    Returns :
        semiaxis_lengths : array of the form [a b] where a is half the length
                           of the axis aligned ellipse along the x-axis and b
                           is half the length along the y-axis (before rotation)
        center : array of the x and y coordinates of the center of the ellipse
        ccw_angle : The counter clockwise angle (in rad) to rotate the ellipse wrt
                    the positive x-axis           
    =#

    S, center = conic2quad(A, B, C, D, E, F);
    return quad2parametric(S, center);
end