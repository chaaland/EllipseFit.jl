module EllipsePlot
    export quadform2ellipse_coords, parametric2ellipse_coords, rotate_mat2d, 
           vandermonde, least_squares_fit, orthogonal_distance_fit, gauss_newton, levenberg_marquardt, conic2parametric

    include("least_squares_fit.jl")
    include("orthogonal_distance_fit.jl")

    ####################################
    # Utils
    ####################################
include("utils/elementwise_pseudoinvert.jl")
    include("utils/utils.jl");    # include the contents of other files in the module
    include("utils/parametric2ellipse_coords.jl");            
    include("utils/quadform2ellipse_coords.jl");
include("utils/ellipse_formatter.jl")
    ####################################
    # Solvers
    #################################### 
    include("solvers/gauss_newton.jl");
    include("solvers/levenberg_marquardt.jl");
    include("solvers/newton_raphson.jl");
    include("solvers/grad_desc.jl");


end
