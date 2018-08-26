
function vandermonde(x::Vector{T}, degree) where T <: Real
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

    N = length(x);
    A = zeros(N, degree + 1);
    col = ones(N);

    for i = 0:degree
        A[:, i+1] = col;
        col = col .* x;
    end

    return A
end
