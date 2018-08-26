# include("../src/utils/quadform2ellipse_coords.jl");
# include("../src/utils/utils.jl");

using PyPlot
using Distributions
using EllipseFit

#= 
Small script demonstrating how to use the module to plot ellipses. 
This script generates a covariance matrix and uses it to plot the various
confidence ellipsoids for a chi-squared distribution with 2 degrees of 
freedom.
=#

pvals = [0.85, 0.9, 0.95, 0.99];

# Generate positive semidefinite covariance matrix
V = rotate_mat2d(pi/3);
D = diagm(0 => vec([3/4 1/4]));
covariance = V * D * V';
precision = inv(covariance);

# Specify mean and degrees of freedom for distribution
mu = [-1 0.5];
dof = 2;

figure(figsize=(10,10))

# Plot confidence ellipses for different false positive rates
for p = pvals
    alpha = cquantile(Chisq(dof), 1-p);         # Probability in right tail
    S = precision / alpha;
    X = quadform2ellipse_coords(S, center=mu);
    plot(X[1,:], X[2,:], label="p = $p");
end

scatter(mu[1], mu[2]);

legend();
title(L"$(x-\mu)^T \Sigma^{-1}(x-\mu) = F_{\chi^2_2}(p)$ confidence ellipsoids");
xlabel(L"$x$");
ylabel(L"$y$");

xlim(-3,3);
ylim(-3,3);
grid(true);

savefig("../img/confidence_ellipsoids.png");
close();