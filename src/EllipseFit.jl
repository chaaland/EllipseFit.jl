module EllipseFit

export ellipse_from_quadform, ellipse_from_parametric, rotate_mat2d, 
        vandermonde, least_squares_fit, orthogonal_distance_fit, gaussnewton, levenbergmarquardt

include("least_squares_fit.jl")
include("orthogonal_distance_fit.jl")

# ####################################
# # Utils
# ####################################
include("utils/utils.jl");    # include the contents of other files in the module
include("utils/ellipse_from_parametric.jl");            
include("utils/ellipse_from_quadform.jl");

# ####################################
# # Solvers
# #################################### 
include("solvers/gaussnewton.jl");
include("solvers/levenbergmarquardt.jl");
include("solvers/newtonraphson.jl");
include("solvers/graddesc.jl");

end