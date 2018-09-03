function [global_B, global_indices, c_OLS_coeff] = create_global_B_matrix(x,...
    y, tri, X, spline_poly_order)
% CREATE_GLOBAL_B_MATRIX Create the sorted global B-matrix for all
% triangles. Also returns the order in which data-points are used and a
% 'symbolic' b-coefficient vector (only used for reference).
% 
% Inputs:
%  - x: x coordinates of meshgrid of vertices
%  - y: y coordinates of meshgrid of vertices
%  - tri: MATLAB triangulation object
%  - X: datapoints
%  - spline_poly_order: order of the polynomial spline
% 
% Outputs:
%  - global_B: global sorted B_form regression matrix
%  - global_indices: order in which datapoints are used in global_B
%  - c_OLS_coeff: b-coefficient vector for later reference
% 
% M.J. Mollema - 03.09.2018

%% Find which data points are in which triangle (not using their barycentric 
% conversion since this did not correspond with the vertex assignment rule)
[IMap, ~] = tsearchn([x(:), y(:)], tri, X);

%% Loop over all triangles to determine B-form regression matrix
% Create empty lists to store data
c_OLS_coeff = [];
global_B = [];
global_indices = [];

% Generate exponentials and coefficients
exponentials    = gen_exp(3, spline_poly_order);
coefficients    = factorial(spline_poly_order) ./ prod(factorial(exponentials), 2);

triangles = sort(tri.ConnectivityList, 2);
for i = 1:size(triangles, 1)
    
    % Get data points and their barycentric coords in current triangle
    indices = find(IMap == i);
    global_indices = vertcat(global_indices, indices);
%     BaryC_current = BaryC(indices, :);
    
    simplex = tri.Points(triangles(i, :), :);
    BaryC_current = bsplinen_cart2bary(simplex, X(indices, :));
    
    % Create sorted B-form regression matrix and global B-form matrix
    sorted_B        = x2fx(BaryC_current, exponentials) .* coefficients';
    global_B        = blkdiag(global_B, sorted_B);
    
    % Create sorted B-coefficient vector for later reference
    c_OLS_coeff           = vertcat(c_OLS_coeff, [i*ones(size(exponentials, 1), 1), exponentials]);
       
end

end