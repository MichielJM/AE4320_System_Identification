%% TODO: create continuity matrix H, h_vector already done but need to add b(v) instead of ones

function do_CompleteSimplex(X, Y, spline_poly_order, spline_cont_order, num_simplices)

order = 0;

% Split data into identification and validation
X_id    = X(2:2:end, 1:2); % Only alpha and beta (???)
X_val   = X(1:2:end, 1:2);
Y_id    = Y(2:2:end);
Y_val   = Y(1:2:end);

% Define triangular grid
grid_start_x    = -0.2;
grid_end_x      = 0.8;
grid_start_y    = -0.3;
grid_end_y      = 0.3;
step_x          = (grid_end_x - grid_start_x)/2^order;
step_y          = (grid_end_y - grid_start_y)/2^order;
[x, y]  = meshgrid(grid_start_x : step_x : grid_end_x,...
                grid_start_y : step_y : grid_end_y);
tri     = delaunayTriangulation(x(:), y(:));
% trimesh(tri, x, y);
% figure;
% plot(x(:), y(:), X(:, 1), X(:, 2))

% Find which data points are in which triangle
[IMap, BaryC] = tsearchn([x(:), y(:)], tri, X_id);
[IMap_val, BaryC_val] = tsearchn([x(:), y(:)], tri, X_val);

% Loop over all triangles to determine B-form regression matrix
c_OLS_coeff = [];
global_B = [];
for i = 1:size(tri.ConnectivityList, 1)
    
    % Get data points and their barycentric coords in current triangle
    indices = find(IMap == i);
    BaryC_current = BaryC(indices, :);
    
    % Create sorted B-form regression matrix and global B-form matrix
    exponentials    = gen_exp(3, spline_poly_order);
    sorted_B        = x2fx(BaryC_current, exponentials);
    global_B        = blkdiag(global_B, sorted_B);
    
    % Create sorted B-coefficient vector for later reference
    c_OLS_coeff           = vertcat(c_OLS_coeff, [i*ones(size(exponentials, 1), 1), exponentials]);
       
end

%% Continuity
% Loop over all orders and edges to determine continuity matrix H
edgess      = edges(tri);
boundary    = freeBoundary(tri);
int_edges   = setdiff(edgess, [boundary; fliplr(boundary)], 'rows'); % Find only interior edges
num_cont_equations = sum(spline_poly_order+1 : -1 : spline_poly_order - spline_cont_order + 1);
H = zeros(num_cont_equations, size(c_OLS_coeff, 1));

for order = 0 : spline_cont_order

    for i = 1 : size(int_edges, 1)
        
        % Left hand part (see lecture 6, slide 76 and onwards for steps)
        exponentials = gen_exp(3, spline_poly_order);
        indices      = exponentials(:, 2) == order;
        LH_part      = exponentials(indices, :);
        
        % Gamma permutations
        gamma = gen_exp(3, order);
        
        % Right hand part
        RH_part = [];
        % Loop over all LH_parts and per part add all gamma permutations
        % and create continuity matrix H
        for j = 1:size(LH_part, 1)
            
            LH_part(:, 2) = 0; % Set second value to 0 (as specified in slide 78)
            LH = LH_part(j, :);
            RH = LH + gamma;
            RH_part = vertcat(RH_part, RH);
            
            % Smoothness conditions (find index of corresponding entry in
            % c_OLS and set that to -1 or b(v))
            h_vector = zeros(1, size(c_OLS_coeff, 1));
            % Left side
            LH(2) = order; % Set
            LH_idx = find(ismember(c_OLS_coeff, [1, LH], 'rows'));
            h_vector(LH_idx) = -1;
            % Right side
            RH_idx = find(ismember(c_OLS_coeff, [2*ones(size(RH, 1), 1), RH], 'rows'));
            h_vector(RH_idx) = 1;
            
            % Fill H
            H(j*(order+1), :) = h_vector;
            
        end        
    end
end

%% Equality strained OLS estimation
Lagrangian = pinv([global_B' * global_B, H';...
            H, zeros(num_cont_equations)]);
C1 = Lagrangian(1:size(H, 2), 1:size(H, 2));
c_OLS = C1 * global_B' * Y_id;

Y_est = c_OLS * X_val;
residuals = Y_est - Y_val;

%% Plotting
% OLS_plotting(X_val, Y_val, Y_est, 1)

end