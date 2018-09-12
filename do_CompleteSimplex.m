function do_CompleteSimplex(X, Y, spline_poly_order, spline_cont_order, ...
    num_triangles_x, num_triangles_y, plot_triangulation, plotting)
% DO_COMPLETE_SIMPLEX Perform data estimation for multiple simplices
% 
% Inputs:
%  - X: Data to estimate with
%  - Y: Data to estimate
%  - spline_poly_order: order of the polynomial spline
%  - spline_cont_order: order to which the splines have to be continuous on
%       triangle edges
%  - num_triangles_x: desired number of triangles in x-direction
%  - num_triangles_y: desired number of triangles in y-direction
%  - plot_triangulation: boolean to determine to plot triangulation
%  - plotting: boolean to determine to plot results 
% 
% Outputs:
%  None
% 
% M.J. Mollema - 09.09.2018

fprintf("\nPerforming complete simplex estimation with order %.f and continuity %.f polynomial\n",...
    spline_poly_order, spline_cont_order)

%% Split data into identification and validation
[X_id, Y_id, X_val, Y_val] = split_data(X, Y, 0.5);

%% Define triangular grid and plot b-coefficient and vertex locations
[tri, x_mesh, y_mesh] = simplex_triangulation(num_triangles_x,...
    num_triangles_y, spline_poly_order, plot_triangulation, X);

%% Create global B matrix
[global_B, global_indices, c_OLS_coeff] = create_global_B_matrix(x_mesh,...
    y_mesh, tri, X_id, spline_poly_order);

%% Continuity
H = simplex_continuity(tri, spline_poly_order, spline_cont_order, c_OLS_coeff);

%% Equality constrained OLS estimation
Lagrangian = pinv([global_B' * global_B, H';...
            H, zeros(size(H, 1), size(H, 1))]);
C1 = Lagrangian(1:size(H,2), 1:size(H,2));
c_OLS = C1 * global_B' * Y_id(global_indices);

% Continuity check
cont_check = H * c_OLS;
if max(cont_check) < 10^(-4)
    fprintf('H*c = 0, continuity assured \n')
else
    fprintf('H*c != 0, continuity not assured \n')
end

%% Validation
% Create global B-matrix for validation
[global_B_val, global_indices_val, ~] = create_global_B_matrix(x_mesh,...
    y_mesh, tri, X_val, spline_poly_order);

% Estimate Y-values
Y_est = global_B_val * c_OLS;

% Variance equal to diagonal of C1
[residual, acx, lags, ~, conf] = model_validation(Y_val(global_indices_val), Y_est, global_B);
VAR = diag(C1);

%% Plotting
if plotting
    % Plot estimation result
    figure; hold on;
    view(20, 12);
    TRI_eval = delaunayn(X_val);
    grey = [192, 192, 192]/255;
    trisurf(TRI_eval, X_val(:, 1), X_val(:, 2), Y_val,...
        'Facecolor', grey, 'FaceAlpha', 0.8, 'EdgeColor', 'None')
    plot3(X_val(global_indices_val, 1), X_val(global_indices_val, 2), Y_est, '.g');
    
    % Plot model residual
    figure;
    plot(residual)
    xlabel("Data point")
    ylabel('Absolute residual')
    
    % Plot model error autocorrelation
    figure; hold on;
    line([lags(1), lags(end)], [conf, conf], 'Color','red','LineStyle','--')
    line([lags(1), lags(end)], [-conf, -conf], 'Color','red','LineStyle','--')
    plot(lags, acx)
    xlabel('Number of lags')
    ylabel('Auto-correlation')
    
    % Plot coefficient variances
    figure; hold on; grid on;
    plot(1:length(VAR), VAR)
    xlabel('Coefficient index')
    ylabel('Coefficient variance')
end

end
