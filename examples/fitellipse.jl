include("noisy_ellipse.jl")
include("../src/least_squares_fit.jl")
include("../src/utils/quadform2ellipse_coords.jl")
include("../src/utils/ellipse_formatter.jl")

using PyPlot

#= 
Small example script demonstrating how to use the module to fit an ellipse
to data. We generate random ellipse data with noise and fit an ellipse to 
the data.
=#

X = noisy_ellipse([3 2], center = [1 -1], ccw_angle=-pi/3, numpoints=100);
A, B, C, D, E, F = least_squares_fit(X);
S, center = conic2quad(A, B, C, D, E, F);

W = quadform2ellipse_coords(S, center=center);

figure(figsize=(10,10))
scatter(X[1,:],X[2,:]);
plot(W[1,:],W[2,:],color="r", label="Least Squares Fit");

xlim([-5,5]);
ylim([-5,5]);
grid(true);

title("Fitting noisy ellipse")
xlabel(L"$x$");
ylabel(L"$y$");
legend();

savefig("../img/fitted_ellipse.png");
close();
