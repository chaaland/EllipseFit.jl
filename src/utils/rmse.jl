
function rmse(x; y=0)
    #= Computes the root mean square error between x and y 

    Given two vectors x and y, return 
        sqrt(1/N * (sum((x_i - y_i)^2)))
    =#

    squared_diff = (x - y).^2;
    mean_square_error = mean(squared_diff);

    return sqrt(mean_square_error)
end