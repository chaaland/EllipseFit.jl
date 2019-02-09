include("solvers.jl")
include("utils/utils.jl")

export Ellipse, QuadraticFormEllipse, ConicFormEllipse, ParametricFormEllipse
export EllipseModel


mutable struct QuadraticFormEllipse{T<:Number}
    S::Array{T,2}
    center::Array{T}

    function QuadraticFormEllipse(S::Array{T,2}, center::Array{T})
        center = vec(center)
        if size(S) != (2,2)
            error("Parameter 'S' must be (2,2)")
        elseif size(center)[1] != 2
            error("Parameter 'center' must be vector of length 2")
        elseif S == zeros(2,2)
            error("Input matrix must be nonzero")
        else
            return new(center, S)
        end
    end

    function QuadraticFormEllipse(S::Array{T,2}) where T<:Number
        return QuadraticFormEllipse(S, vec([0 0]))
    end
end
QuadraticFormEllipse(S::Array{Real,2}, center::Array{Real}) = Ellipse(promote(S, center))

mutable struct ConicFormEllipse{T<:Number}
    A::T
    B::T
    C::T
    D::T
    E::T

    function ConicFormEllipse(A::T, B::T, C::T, D::T, E::T)
        if B ^2 - 4 * A * C >= 0
            error("Discriminant is non-negative. Input does not denote an ellipse")
        else
            return new(A, B, C, D, E)
        end
    end
end
ConicFormEllipse(A::Real, B::Real, C::Real, D::Real, E::Real) = ConicFormEllipse(promote(A, B, C, D, E))

mutable struct ParametricFormEllipse{T<:Real}
    center::Array{T}
    semiaxis_lengths::Array{T}
    ccw_angle::T

    function ParametricFormEllipse(center::Array{T}, semiaxis_lengths::Array{T}, ccw_angle::Array{T})
        center = vec(center)
        semiaxis_lengths = vec(semiaxis_lengths)

        if size(center)[1] != 2
            error("Parameter 'center' must be vector of length 2")
        elseif size(semiaxis_lengths)[1] != 2
            error("Parameter 'semiaxis_lengths' must be vector of length 2")
        elseif sum(semiaxis_lengths .< 0) > 0
            error("Parameter 'semiaxis_lengths' must be nonnegative")
        else
            return new(center, semiaxis_lengths, ccw_angle)
        end
    end
end
ParametricFormEllipse(center::Array{Real}, semiaxis_lengths::Array{Real}, ccw_angle) = ParametricFormEllipse(promote(center, semiaxis_lengths), ccw_angle)

mutable struct Ellipse
    quadratic::QuadraticFormEllipse
    conic::ConicFormEllipse
    parametric::ParametricFormEllipse
end

function Ellipse(center::Array{T}, S::Array{T,2}) where T<:Real
    quadform = QuadraticFormEllipse(center, S)
    conicform = quad2conic(quadform)
    parametricform = quad2parametric(quadform)
    return Ellipse(quadform, conicform, parametricform)
end
Ellipse(center::Array{Real}, S::Array{Real,2}) = Ellipse(promote(center, S))

function Ellipse(A::T, B::T, C::T, D::T, E::T) where T<: Real
    conicform = ConicFormEllipse(A, B, C, D, E)
    quadform = conic2quad(conicform)
    parametricform = conic2parametric(conicform)
    return Ellipse(quadform, conicform, parametricform)
end
Ellipse(A::Real, B::Real, C::Real, D::Real, E::Real) = Ellipse(promote(A, B, C, D, E))

function Ellipse(center::Array{T}, semiaxis_lengths::Array{T}, ccw_angle=0) where T<:Real
    center, semiaxis_lengths = promote(center, semiaxis_lengths )
    parametricform = ParametricFormEllipse(center, semiaxis_lengths, ccw_angle)
    quadform = parametric2quad(quadform)
    conicform = parametric2conic(parametricform)
    return Ellipse(quadform, conicform, parametricform)
end
Ellipse(center::Array{Real}, semiaxis_lengths::Array{Real}, ccw_angle=0) = Ellipse(promote(center, semiaxis_lengths), ccw_angle)

struct EllipseModel
    X::Array{Real,2}
    objective::Objective
    solver::Solver 
    solution::Union{Ellipse,Nothing}

    function EllipseModel(X::Array{T,2}, objective::Objective, solver::Solver, 
                           solution::Union{Ellipse,Nothing}=nothing) where T <: Real
        N, n = size(X)
        if n != 2
            error("Expected an array with second dimension 2")
        end
        new(X, objective, solver, solution)
    end
end

function quad2conic(qform::QuadraticFormEllipse)
   #= Helper for converting quadratic form ellipse into a standard conic form
    
    Given an ellipse as a quadratic form

            (x - c)^T S (x - c) = 1

    convert it to the conic section form 

        A * x^2 + B * x * y + C * y^2 + D * x + E * y + F = 0

    Args :
        qform : quadratic form ellipse to be converted to conic form

    Returns :
        ellipse in conic form
        
    =# 
    S = qform.S
    center = qform.center

    F = center' * S * center - 1;

    A = S[1,1]
    B = (S[1,2] + S[2,1]) / F
    C = S[2,2]

    b = -2 * center' * S;
    D = b[1]
    E = b[2]

    if (F == 0)
        ConicFormEllipse(A, B, C, D, E)
    end

    return ConicFormEllipse(A / F, B / F, C / F, D / F, E / F)
end

function quad2parametric(qform::QuadraticFormEllipse)
   #= Helper for converting from quadratic form to parameteric form of ellipse

    Given an ellipse as a quadratic form 

            (x - c)^T S (x - c) = 1

    convert it to the parametric form

            [x_c y_c] + rotation_mat(ccw_angle) * [a*cos(theta) b*sin(theta)]

    Args :
        qform : quadratic form ellipse to be converted to conic form 

    Returns :
        parametric form ellipse
    =#

    f = eigen(qform.S);
    V = f.vectors;
    D = f.values;

    semiaxis_lengths = sqrt(elementwise_pseudoinvert(D));
    ccw_angle = atan(V[2,1], V[1,1]);

    return ParametricFormEllipse(qform.center, semiaxis_lengths, ccw_angle)
end

function conic2quad(cform::ConicFormEllipse)
    #= Helper for converting conic form ellipse into a standard quadratic form
    
    Given an ellipse as a conic section of the form

        A * x^2 + B * x * y + C * y^2 + D * x + E * y + F = 0

    convert it to the form 

        (x - c)^T S (x - c) = 1

    This is done by completing the square leading to the following equivalence
            x^T * Q * x + b^T * x + c = 0
        =>  (x + 0.5 * inv(Q) * b)^T Q (x + 0.5 * inv(Q) * b)) + c - 0.25 * b^T * inv(A) * b = 0
    
    Args :

    Returns :

    =#
    A = cform.A
    B = cform.B
    C = cform.C
    D = cform.D
    E = cform.E

    Q = [A B/2; B/2 C];
    b = [D; E];

    beta = Q \ b;
    c = 0.25 * b' * beta;

    S = Q / c;
    center = -0.5 * beta;

    return QuadraticFormEllipse(S, center)
end

function conic2parametric(cform::ConicFormEllipse)
    #= Helper for converting from conic form to parametric form of ellipse

    Given an ellipse in standard conic section form

            A * x^2 + B * x * y + C * y^2 + D * x + E * y + F = 0

    convert it to parametric form

            [x_c y_c] + rot_mat2d(ccw_angle) * [a*cos(theta) b*sin(theta)]

    Args :
        

    Returns :
        semiaxis_lengths : array of the form [a b] where a is half the length
                           of the axis aligned ellipse along the x-axis and b
                           is half the length along the y-axis (before rotation)
        center : array of the x and y coordinates of the center of the ellipse
        ccw_angle : The counter clockwise angle (in rad) to rotate the ellipse wrt
                    the positive x-axis           
    =#

    qform = conic2quad(cform);
    return quad2parametric(qform)
end

function parametric2quad(pform::ParametricFormEllipse)
    #= Helper for converting from paramteric to quadratic form of ellipse

    Given an ellipse in parametric form

        [x_c y_c] + rotation_mat(ccw_angle) * [a*cos(theta) b*sin(theta)]

    covnert it into a quadratic of the form

        (x - c)^T S (x - c) = 1

    Args :
        
    Returns :
       
    =#

    sqrtD = diagm(0 => vec(pform.semiaxis_lengths));
    invD = elementwise_pseudoinvert(sqrtD .^ 2);
    V = rotation_mat(pform.ccw_angle);
    S = V * invD * V';

    return QuadraticFormEllipse(S, pform.center)
end

function parametric2conic(pform::ParametricFormEllipse)
    #= Helper for converting from parametric form to standard form of ellipse

    Given an ellipse in parametric form

            [x_c y_c] + rotation_mat(ccw_angle) * [a*cos(theta) b*sin(theta)]

    convert it to the standard conic section form

            A * x^2 + B * x * y + C * y^2 + D * x + E * y + F = 0

    Args :
        

    Returns :
              
    =#

    qform = parametric2quad(pform)
    return conic2quad(qform)
end