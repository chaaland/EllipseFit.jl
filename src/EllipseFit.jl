module EllipseFit

include("solvers.jl")
include("ellipse.jl")
# include("fit.jl")

####################################
# Utils
####################################
include("utils.jl");    # include the contents of other files in the module

####################################
# Solvers
#################################### 
include("solvers/gauss_newton.jl");
include("solvers/levenberg_marquardt.jl");
include("solvers/newton_raphson.jl");
include("solvers/grad_desc.jl");

end