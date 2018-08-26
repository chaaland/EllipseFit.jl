
function orthogonal_dist_fit(X::Array{T,2}) where T <: Real
    #= Fit an ellipse by minimizing the orthogonal distance of the points 

    Rather than least squares, the ellipse is fit so as to minimize the 
    perpendicular distance of all the measured points an ellipse. This is
    an example of an errors in variables model which accounts for measurement
    error in both the independent and dependent variables

    Args :
        X : 2 x N or N x 2 matrix containing the data to be fit with an ellipse

    Returns :

    =#
    m, n = size(X);

    if m == 2  
        X = 1.0 * X';
    end
    N = size(X, 1);
    num_params = 5 + N;             # xcenter, ycenter, angle, semi major, semi minor + theta per data point
    num_equalities = 2 * N;
    thetavals_lm, fvals_lm, gradnorm_lm, lambdavals_lm = levenberg_marquardt((num_params, num_equalities), parametric_ellipse, 
                                                                    jacobian_ellipse; xinit=xinit, max_iters=100, atol=1e-6);

    center = vec(thetavals_lm[1:2,end]);
    semiaxis_lengths = vec(thetavals_lm[3:4,end]);
    ccw_angle = thetavals_lm[5, end];

    return center, semiaxis_lengths, ccw_angle
end

function parametric_ellipse(x::Array{T,1}) where T <: Real
    center = x[1:2];
    semiaxis_lengths = x[3:4];
    ccw_angle = x[5];
    theta = x[6:end];
    
    onaxis_ellipse = Diagonal(semiaxis_lengths) * [cos.(theta) sin.(theta)]';
    rotated_ellipse = rotate_mat2d(ccw_angle) * onaxis_ellipse;
    shifted_ellipse = center .+ rotated_ellipse;

    return vec(shifted_ellipse) - vec(X')
end

function jacobian_ellipse(x::Array{T,1}) where T <: Real
    center = vec(x[1:2]);
    semiaxis_lengths = vec(x[3:4]);
    ccw_angle = x[5];
    theta = vec(x[6:end]);
    
    a = semiaxis_lengths[1];
    b = semiaxis_lengths[2];
    
    dxc = repeat([1, 0], length(theta), 1);
    dyc = repeat([0, 1], length(theta), 1);
    da = vec([cos.(ccw_angle)*cos.(theta) sin.(ccw_angle)*cos.(theta)]');
    db = vec([-sin.(ccw_angle)*sin.(theta) cos.(ccw_angle)*sin.(theta)]');
    dalpha = vec([-a*sin.(ccw_angle)*cos.(theta)-b*cos.(ccw_angle)*sin.(theta) a*cos.(ccw_angle)*cos.(theta)-b*sin.(ccw_angle)*sin.(theta)]');
    dtheta = zeros(2*length(theta), length(theta));
    
    for i = 1:length(theta) 
        theta_val = theta[i];
        val1 = -a * cos.(ccw_angle) * sin.(theta_val) - b * sin.(ccw_angle) * cos.(theta_val);
        val2 = -a * sin.(ccw_angle) * sin.(theta_val) + b * cos.(ccw_angle) * cos.(theta_val);
        dtheta[2*i-1:2*i, i] = [val1 val2];
    end

    J = [dxc dyc da db dalpha dtheta];
end
