module EllipsePlot
    export quadform2ellipse_coords, parametric2ellipse_coords, rotate_mat2d, 
           vandermonde, gauss_newton, levenberg_marquardt

    ####################################
    # Utils
    ####################################
    include("utils/elementwise_pseudoinvert.jl");    # include the contents of other files in the module
    include("utils/rotate_mat2d.jl");
    include("utils/parametric2ellipse_coords.jl");            
    include("utils/quadform2ellipse_coords.jl");
    include("utils/vandermonde.jl");

    ####################################
    # Solvers
    #################################### 
    include("solvers/gauss_newton.jl");
    include("solvers/levenberg_marquardt.jl");
    include("solvers/newton_raphson.jl");
    include("solvers/grad_desc.jl");

end