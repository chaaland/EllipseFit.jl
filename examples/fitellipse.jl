using PyPlot
using Random
include("../src/EllipseFit.jl")
using .EllipseFit

#=
Small example script demonstrating how to use the module to fit an ellipse
to data. We generate random ellipse data with noise and fit an ellipse using both least
squares and orthogonal distance regression
=#

## generate noisy ellipse
Random.seed!(3084021)
n_plot_points = 500
X = noisyellipse([3 2], center = [1 -1], ccw_angle=-pi/3, numpoints=100);

## Fit
xinit = vcat([0, 0, 1, 1, 0], randn(size(X,1)));
ellipse_model_ls = EllipseModel(X, LeastSquares, NormalEquations())
ellipse_model_nlls = EllipseModel(X, OrthogonalEuclideanDistance, LevenbergMarquardt(xinit=xinit))
ls_fit_ellipse = fit(ellipse_model_ls)
nlls_fit_ellipse = fit(ellipse_model_nlls)
X_fitted_ls = ellipse_to_plot_points(ls_fit_ellipse) #, n=n_plot_points)
X_fitted_nlls = ellipse_to_plot_points(nlls_fit_ellipse) #, n=n_plot_points)

## Largests residuals
n_res = 10;
ls_res = get_residual(ls_fit_ellipse, X);
ind = partialsortperm(abs.(ls_res),1:n_res,rev=true);
largest_res_ls = X[ind,:];
nlls_res = get_residual(nlls_fit_ellipse, X);
ind = partialsortperm(abs.(nlls_res),1:n_res,rev=true);
largest_res_nlls = X[ind,:];

## Plot results
fig = figure(figsize=(7,7))

scatter(X[:,1], X[:,2]);
scatter(largest_res_ls[:,1],largest_res_ls[:,2],marker=:x,label="LS largest errors")
scatter(largest_res_nlls[:,1],largest_res_nlls[:,2],marker=:+, label="NLLS largest errors")
plot(X_fitted_ls[:,1], X_fitted_ls[:,2], color="r", label="Least Squares Fit");
plot(X_fitted_nlls[:,1], X_fitted_nlls[:,2], color="m", label="Nonlinear Least Squares Fit");

xlim([-5,5]);
ylim([-5,5]);
grid(true);

title("Fitting noisy ellipse")
xlabel(L"$x$");
ylabel(L"$y$");
legend();
display(fig)
# savefig("img/fitted_ellipse.png");
