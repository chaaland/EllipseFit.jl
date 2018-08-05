
function elementwise_pseudoinvert(v::AbstractArray, tol=1e-10) 
    #= Helper for elementwise inversion of a matrix v
    
    Inverts the non-zero elements of v and keeps the 0 elements. 
    Elements within tol of 0 are treated as zero
    
    Args :
        v : array
        tol : 
    
    Returns :
        Array with the elements inverted except the zeros
    =#

    m = maximum(abs.(v));
    v = v ./ m;
    reciprocal = 1./v;
    reciprocal[abs.(reciprocal) .>= 1/tol] = 0;   

    return reciprocal / m;
end
