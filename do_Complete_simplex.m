function do_Complete_simplex_REDO(X, Y, spline_poly_order, spline_cont_order, ...
    num_triangles_x, num_triangles_y, plot_triangulation)
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
% 
% Outputs:
%  None

%% Split data into identification and validation
% Even indices for identification, odd indices for validation
X_id    = X(2:2:end, 1:2);
X_val   = X(1:2:end, 1:2);
Y_id    = Y(2:2:end);
Y_val   = Y(1:2:end);

%% Define triangular grid and plot b-coefficient and vertex locations
[tri, x_mesh, y_mesh] = simplex_triangulation(num_triangles_x,...
    num_triangles_y, spline_poly_order, plot_triangulation);

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

%% Validation
% Create global B-matrix for validation
[global_B_val, global_indices_val, ~] = create_global_B_matrix(x_mesh,...
    y_mesh, tri, X_val, spline_poly_order);

% Estimate Y-values
Y_est = global_B_val * c_OLS;

%% Plotting  
figure; hold on;
TRI_eval = delaunayn(X_val);
grey = [192, 192, 192]/255;
trisurf(TRI_eval, X_val(:, 1), X_val(:, 2), Y_val,...
    'Facecolor', grey, 'FaceAlpha', 0.8, 'EdgeColor', 'None')
plot3(X_val(global_indices_val, 1), X_val(global_indices_val, 2), Y_est, '.g');
% plot3(X_val(:, 1), X_val(:, 2), Y_val, '.k');

end