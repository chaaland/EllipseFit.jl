using LinearAlgebra; 

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

    semiaxis_lengths = sqrt.(elementwise_pseudoinvert(D));
    ccw_angle = acos(V[1,1]);

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

function conic2area(A, B, C, D, E, F)
    #= Return area of an ellipse in standard form
    Given an ellipse in standard conic section form

            A * x^2 + B * x * y + C * y^2 + D * x + E * y + F = 0

    Return its area
    Taken from https://math.stackexchange.com/questions/247332/area-of-an-ellipse

    When rotating conics in implicit form
    $$
    Ax^2+Bxy+Cy^2+Dx+Ey+F=0\tag{1}
    $$
    around the origin there are 5 invariants:
    $$
    \begin{array}{rl}
    I_1&=A+C\\
    I_2&=(A-C)^2+B^2\\
    I_3&=D^2+E^2\\
    I_4&=(A-C)(D^2-E^2)+2DEB\\
    I_5&=F\tag{2}
    \end{array}
    $$
    Assuming that we have rotated to eliminate $B$, we have
    $$
    Ax^2+Cy^2+Dx+Ey+F=0\tag{3}
    $$
    Translating the center to the origin gives
    $$
    Ax^2+Cy^2+F-\frac{D^2}{4A}-\frac{E^2}{4C}=0\tag{4}
    $$
    which is the same as
    $$
    \frac{x^2}{a^2}+\frac{y^2}{b^2}=1\tag{5}
    $$
    if we set
    $$
    a=\sqrt{\frac{\frac{D^2}{4A}+\frac{E^2}{4C}-F}{A}}\quad\text{and}\quad
    b=\sqrt{\frac{\frac{D^2}{4A}+\frac{E^2}{4C}-F}{C}}\tag{6}
    $$

    Using $\text{Area}=\pi ab$ and $(4)$ and rewriting in terms of the invariants to remove the rotation, we get
    $$
    \begin{align}
    \text{Area}
    &=\pi\frac{\frac{D^2}{4A}+\frac{E^2}{4C}-F}{\sqrt{AC}}\\
    &=\pi\frac{\color{#C00000}{2CD^2+2AE^2}-\color{#00A000}{8ACF}}{\color{#0000FF}{8(AC)^{3/2}}}\\
    &=\pi\frac{\color{#C00000}{I_1I_3-I_4}-\color{#00A000}{2(I_1^2-I_2)I_5}}{\color{#0000FF}{(I_1^2-I_2)^{3/2}}}\tag{7}
    \end{align}
    $$
    Therefore, $(2)$ and $(7)$ give the area in terms of the coefficients in $(1)$.

    **Invariants Under Rotation and Translation**

    $I_1$ and $I_2$ are invariant under rotation and translation, but there is one more: the constant coefficient when the center is translated to the origin. Writing $F-\dfrac{D^2}{4A}-\dfrac{E^2}{4C}$ in terms of the rotational invariants and expanding yields
    $$
    I_6=F-\frac{AE^2-BDE+CD^2}{4AC-B^2}\tag{8}
    $$
    Thus, the conic satisfying $(1)$, when rotated and translated becomes
    $$
    (I_1-\sqrt{I_2})x^2+(I_1+\sqrt{I_2})y^2+2I_6=0\tag{9}
    $$

    Args :
        A : Coefficient of x^2
        B : Coefficient of x * y 
        C : Coefficient of y^2
        D : Coefficient of x
        E : Coefficient of y 
        F : Bias 

    Returns : area
    =#
    return pi * (2*C*D^2 + 2*A*E^2 - 8 * A * C * F)/(8 * (A*C)^(3/2))
end
