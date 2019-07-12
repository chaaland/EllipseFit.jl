using PyPlot
using EllipseFit

# include("../utils.jl")
#= 
Small example script demonstrating how to use the module to fit an ellipse
to data. We generate random ellipse data with noise and fit an ellipse to 
the data.
=#

n_plot_points = 500
X = noisyellipse([3 2], center = [1 -1], ccw_angle=-pi/3, numpoints=100);
ellipse_model_ls = EllipseModel(X, LeastSquares, NormalEquations)
ellipse_model_nnls = EllipseModel(X, OrthogonalEuclideanDistance, LevenbergMarquardt)
ls_fit_ellipse = fit(ellipse_model_ls)
nnls_fit_ellipse = fit(ellipse_model_nnls)
X_fitted_ls = ellipse_to_plot_points(ls_fit_ellipse) #, n=n_plot_points)
X_fitted_nnls = ellipse_to_plot_points(nnls_fit_ellipse) #, n=n_plot_points)

fig = figure(figsize=(7,7))
print(X[:,1])
scatter(X[:,1], X[:,2]);
plot(X_fitted_ls[:,1], X_fitted_ls[:,2],color="r", label="Least Squares Fit");
plot(X_fitted_nnls[:,1], X_fitted_nnls[:,2],color="r", label="Nonlinear Least Squares Fit");

xlim([-5,5]);
ylim([-5,5]);
grid(true);

title("Fitting noisy ellipse")
xlabel(L"$x$");
ylabel(L"$y$");
display(fig)
savefig("img/fitted_ellipse.png");
legend();
