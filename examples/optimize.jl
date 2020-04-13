
using PyPlot
using Random
using EllipseFit

Random.seed!(1)

z = vec(rand(-5:5,2));
A = vandermonde(randn(10), 1);
b = A * z + vec(randn(10));


function f(x)
    r = A * x - b;
    return sum(r.^2);
end

function grad(x)
    return 2 * A' * (A * x - b);
end

function jacobian(x)
    return grad(x)'
end

function hessian(x)
    return 2 * A' * A;
end

traj_gd, cost_gd, gradmagnitude_gd = graddesc(2, f, grad, max_iters=100, alpha=0.01);
traj_nr, cost_nr, gradmagnitude_nr = newtonraphson(2, f, jacobian, hessian, max_iters=10);

figure(figsize=(10,10));
x = vec(LinRange(-5,3,100));
y = vec(LinRange(-3,8,100));
    
X = repmat(x, 1, 100);
Y = repmat(y', 100, 1);

Z = [norm(A * vec([i j]) - b)^2 for i in x, j in y];
contour(X, Y, Z, 15);
plot(traj_gd[1,:], traj_gd[2,:], marker="o", label="GradientDescent");
plot(traj_nr[1,:], traj_nr[2,:], marker="o", label="NewtonRaphson");
scatter(traj_gd[1,end], traj_gd[2,end], marker="*");
legend();
savefig("img/optimize.png")