module EllipsePlot
    export quadform2ellipse_coords, parametric2ellipse_coords 

    ####################################
    # Utils
    ####################################
    include("utils/elementwise_pseudoinvert.jl");    # include the contents of other files in the module
    include("utils/rotate_mat2d.jl");
    include("utils/parametric2ellipse_coords.jl");            
    include("utils/quadform2ellipse_coords.jl");

end