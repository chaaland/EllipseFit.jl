# Small script demonstrating how to use the module to plot ellipses. 
# Generate a covariance matrix and uses it to plot confidence ellipsoids 
# for a chi-squared distribution with 2 degrees of freedom.

using PyPlot
using Distributions
using Random
using EllipseFit

# Set seed for reproducibility
Random.seed!(3084021)

# Set some false positive rates (1 - p)
pvals = [0.85, 0.9, 0.95, 0.99];

# Generate positive semidefinite covariance matrix
dof = 2;
R = rand(dof, dof);
covariance = R' * R;
precision = inv(covariance);

# Specify mean and degrees of freedom for distribution
mu = [-1 0.5];

fig = figure(figsize=(7,7));

# Plot confidence ellipses 
for p = pvals
    alpha = cquantile(Chisq(dof), (1 - p)/2.0);         # Probability in right tail
    S = precision / alpha;
    
    confidence_ellipse = Ellipse(S, center=mu)
    X_plot = ellipse_to_plot_points(confidence_ellipse) #, n=n_plot_points)
    plot(X_plot[:,1], X_plot[:,2], label="p = $p");
end
scatter(mu[1], mu[2]);

legend();
title(L"$(x-\mu)^T \Sigma^{-1}(x-\mu) = F_{\chi^2_2}(p)$ confidence ellipsoids",fontsize=14);
xlabel(L"$x$",fontsize=14);
ylabel(L"$y$",fontsize=14);

grid(true, which="major");
grid(true, which="minor",linestyle="--");
PyPlot.minorticks_on();
display(fig)
savefig("img/confidence_ellipsoids.png");