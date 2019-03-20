using LinearAlgebra
include("utils.jl")

export Ellipse, QuadraticFormEllipse, ConicFormEllipse, ParametricFormEllipse

struct QuadraticFormEllipse{T<:Real, U<:Real}
    S::Array{T,2}
    center::Array{U}

    function QuadraticFormEllipse(S::Array{T,2}, center::Array{U}) where {T<:Real,U<:Real}
        center = vec(center)
        if size(S) != (2,2)
            error("Parameter 'S' must be (2,2)")
        elseif size(center)[1] != 2
            error("Parameter 'center' must be vector of length 2")
        elseif S == zeros(2,2)
            error("Input matrix must be nonzero")
        end
        new{T,U}(S, center)
    end
    QuadraticFormEllipse(x::Array{T,2}) where {T<:Real} = QuadraticFormEllipse(x, [0.0 0.0])
end

struct ConicFormEllipse{T<:Real}
    A::T
    B::T
    C::T
    D::T
    E::T

    function ConicFormEllipse(A::T, B::T, C::T, D::T, E::T) where {T<:Real}
        if B ^2 - 4 * A * C >= 0
            error("Discriminant is non-negative. Input does not denote an ellipse")
        end

        new{T}(A, B, C, D, E)
    end
end
ConicFormEllipse(A::Real, B::Real, C::Real, D::Real, E::Real) = ConicFormEllipse(promote(A, B, C, D, E)...)

struct ParametricFormEllipse{T<:Real, U<:Real, V<:Real}
    semiaxis_lengths::Array{T}
    center::Array{U}
    ccw_angle::V

    function ParametricFormEllipse(semiaxis_lengths::Array{T}, center::Array{U}, ccw_angle::V) where {T<:Real, U<:Real, V<:Real}
        semiaxis_lengths = vec(semiaxis_lengths)
        center = vec(center)

        if size(center)[1] != 2
            error("Parameter 'center' must be vector of length 2")
        elseif size(semiaxis_lengths)[1] != 2
            error("Parameter 'semiaxis_lengths' must be vector of length 2")
        elseif sum(semiaxis_lengths .< 0) > 0
            error("Parameter 'semiaxis_lengths' must be nonnegative")
        else
            new{T,U,V}(semiaxis_lengths, center, ccw_angle)
        end
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
 
    A = S[1,1]
    B = (S[1,2] + S[2,1])
    C = S[2,2]
 
    b = -2 * center' * S;
    D = b[1]
    E = b[2]
    negF = 1 - center' * S * center
 
    return ConicFormEllipse(A/negF, B/negF, C/negF, D/negF, E/negF)
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
 
    f = eigen(qform.S)
    V = f.vectors
    D = f.values                
 
    semiaxis_lengths = sqrt.(elementwise_pseudoinvert(D))
    p = sortperm(semiaxis_lengths, rev=true)
    sorted_semiaxes = semiaxis_lengths[p]
    sorted_eig_vecs = V[:,p]
    major_axis = sorted_eig_vecs[:,1]
    ccw_angle = atan(major_axis[2], major_axis[1])
 
    return ParametricFormEllipse(sorted_semiaxes, qform.center, ccw_angle)
end

function conic2quad(cform::ConicFormEllipse)
    #= Helper for converting conic form ellipse into a standard quadratic form
    
    Given an ellipse as a conic section of the form

        A * x^2 + B * x * y + C * y^2 + D * x + E * y + F = 0

    convert it to the form 

        (x - c)^T S (x - c) = 1

    This is done by completing the square leading to the following equivalence
            x^T * Q * x + b^T * x + c = 0
        =>  (x + 0.5 * inv(Q) * b)^T Q (x + 0.5 * inv(Q) * b)) + c - 0.25 * b^T * inv(Q) * b = 0
    
    Args :

    Returns :

    =#

    A = cform.A
    B = cform.B
    C = cform.C
    D = cform.D
    E = cform.E

    Q = [A B/2; B/2 C]
    b = [D; E]

    beta = Q \ b
    rhs = 1 + 0.25 * b' * beta

    S = Q / rhs
    center = -0.5 * beta

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

    qform = conic2quad(cform)
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

    sqrtD = diagm(0 => vec(pform.semiaxis_lengths))
    invD = elementwise_pseudoinvert(sqrtD .^ 2)
    V = rotation_mat(pform.ccw_angle)
    S = V * invD * V'

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
    return quad2conic(qform)
end 

struct Ellipse
    quadform::QuadraticFormEllipse
    conicform::ConicFormEllipse
    parametricform::ParametricFormEllipse
end

function Ellipse(S::Array{T,2}, center::Array{U}) where {T<:Real, U<:Real}
    quadform = QuadraticFormEllipse(S, center)
    conicform = quad2conic(quadform)
    parametricform = quad2parametric(quadform)
    return Ellipse(quadform, conicform, parametricform)
end

function Ellipse(A::T, B::T, C::T, D::T, E::T) where T<: Real
    conicform = ConicFormEllipse(A, B, C, D, E)
    quadform = conic2quad(conicform)
    parametricform = conic2parametric(conicform)
    return Ellipse(quadform, conicform, parametricform)
end
Ellipse(A::Real, B::Real, C::Real, D::Real, E::Real) = Ellipse(promote(A, B, C, D, E)...)

function Ellipse(semiaxis_lengths::Array{T}; center=[0 0], ccw_angle=0) where T<:Real
    parametricform = ParametricFormEllipse(semiaxis_lengths, center, ccw_angle)
    quadform = parametric2quad(parametricform)
    conicform = parametric2conic(parametricform)
    return Ellipse(quadform, conicform, parametricform)
end
