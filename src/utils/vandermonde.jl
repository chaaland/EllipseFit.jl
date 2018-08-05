
function vandermonde(x::Array{T,1}, n) where T <: Real
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

    A = zeros(length(x), n + 1);
    xtilde = vec(x);
    col = ones(size(xtilde));
    for i = 0:n
        A[:,i+1] = col;
        col = col .* xtilde;
    end

    return A
end