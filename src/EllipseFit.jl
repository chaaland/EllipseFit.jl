module EllipseFit

include("solvers.jl")
include("ellipse.jl")
# include("fit.jl")

####################################
# Utils
####################################
include("utils.jl")    # include the contents of other files in the module

####################################
# Solvers
#################################### 
include("solvers/gaussnewton.jl")
include("solvers/levenbergmarquardt.jl")
include("solvers/newtonraphson.jl")
include("solvers/graddesc.jl")

end