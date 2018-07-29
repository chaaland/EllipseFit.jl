include("../src/utils/quadform2ellipse_coords.jl");
include("../src/utils/rotate_mat2d.jl");

using PyPlot
using Distributions

#= 
Small script demonstrating how to use the module to plot ellipses. 
This script generates a covariance matrix and uses it to plot the various
confidence ellipsoids for a chi-squared distribution with 2 degrees of 
freedom.
=#

pvals = [0.85, 0.9, 0.95, 0.99];

V = rotate_mat2d(pi/3);
D = diagm(vec([3/4 1/4]));
covariance = V * D * V';
precision = inv(covariance);
print(precision ≈ precision')
mu = [-1 0.5];
dof = 2;

figure(figsize=(10,10))
for p = pvals
    alpha = cquantile(Chisq(dof), 1-p);         # Probability in right tail
    S = precision / alpha;
    X = quadform2ellipse_coords(S, center=mu);
    plot(X[1,:], X[2,:], label="p = $p");
end

scatter(c[1], c[2]);

title(L"$(x-\mu)^T \Sigma^{-1}(x-\mu) = F_{\chi^2_2}(p)$ confidence ellipsoids");
xlabel(L"$x$");
ylabel(L"$y$");

legend();
xlim(-3,3);
ylim(-3,3);
grid(true);

savefig("../img/confidence_ellipsoids.png");
close();