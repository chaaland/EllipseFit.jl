module EllipseFit

include("solvers.jl")
include("ellipse.jl")
include("fit.jl")

####################################
# Utils
####################################
include("utils.jl")

####################################
# Solvers
#################################### 
include("solvers.jl")
include("solvers/gaussnewton.jl")
include("solvers/levenbergmarquardt.jl")
include("solvers/newtonraphson.jl")
include("solvers/graddesc.jl")

end