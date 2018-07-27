module EllipsePlot
    export quadform2ellipse_coords, parametric2ellipse_coords 

    ####################################
    # Utils
    ####################################
    include("utils/elementwise_pseudoinvert.jl");    # include the contents of other files in the module
    include("utils/rotate_mat2d.jl");
    include("parametric2ellipse_coords.jl");            
    include("quadform2ellipse_coords.jl");

end