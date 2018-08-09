%% TODO: create continuity matrix H, h_vector already done but need to add b(v) instead of ones

function do_CompleteSimplex(X, Y, spline_poly_order, spline_cont_order, num_simplices)

order = 1;

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
trimesh(tri, x, y);
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
% H = zeros(num_cont_equations * size(int_edges, 1), size(c_OLS_coeff, 1));
H = [];

for order = 0 : spline_cont_order
    order
    for i = 1 : size(int_edges, 1)
        
        % Find triangles and their out-of-edge vertices connected to current edge
        triangle_IDs    = edgeAttachments(tri, int_edges(i, :));
        vertex_IDs      = tri.ConnectivityList(triangle_IDs{:}, :);
        vertex_ID_1     = setdiff(vertex_IDs(1, :), vertex_IDs(2, :)); % Find values in A which are not in B
        vertex_ID_2     = setdiff(vertex_IDs(2, :), vertex_IDs(1, :)); % Find values in B which are not in A
        vertex_cart_1   = tri.Points(vertex_ID_1, :);
        vertex_cart_2   = tri.Points(vertex_ID_2, :);
        % Determine barycentric coords of out-of-edge vertices wrt their
        % own triangle and the other triangle (11 is own, 12 wrt other etc)
        vertex_bary_11  = cartesianToBarycentric(tri, triangle_IDs{1}(1), vertex_cart_1);
        vertex_bary_12  = cartesianToBarycentric(tri, triangle_IDs{1}(2), vertex_cart_1);
        vertex_bary_21  = cartesianToBarycentric(tri, triangle_IDs{1}(1), vertex_cart_2);
        vertex_bary_22  = cartesianToBarycentric(tri, triangle_IDs{1}(2), vertex_cart_2);        
        
        % Left hand part (see lecture 6, slide 76 and onwards for steps)
        single_nonzero_loc_1 = find(vertex_bary_11 == 1);
        single_nonzero_loc_2 = find(vertex_bary_22 == 1);
        exponentials = gen_exp(3, spline_poly_order);
        indices      = exponentials(:, single_nonzero_loc_1) == order;  %TODO: MAKE THIS VALUE DEPENDENT ON THE SINGLE NONZERO VALUE OF OUT OF EDGE VERTEX
        LH_part      = exponentials(indices, :);
        
        % Gamma permutations
        gamma = gen_exp(3, order);
        
        % Right hand part
        RH_part = [];
        % Loop over all LH_parts and per part add all gamma permutations
        % and create continuity matrix H
        for j = 1:size(LH_part, 1)
            
            % Get k0 and k1 from left hand part
            ks = LH_part(j, :);
            ks(:, single_nonzero_loc_1) = [];
            RH = zeros(1, 3);
            if single_nonzero_loc_2 == 1
                RH(2:3) = ks;
            elseif single_nonzero_loc_2 == 2
                RH(1) = ks(1);
                RH(3) = ks(2);
            else
                RH(1:2) = ks;
            end
            
            % All permutations of RH + gamma
            RH = RH + gamma;
            RH_part = vertcat(RH_part, RH);
            
%             LH_part(:, single_nonzero_loc_2) = 0 % slide 78 TODO: MAKE THIS VALUE DEPENDENT ON THE SINGLE NONZERO VALUE OF OUT OF EDGE VERTEX
%             LH = LH_part(j, :);
%             RH = LH + gamma;
%             RH_part = vertcat(RH_part, RH);
%             
            % Smoothness conditions (find index of corresponding entry in
            % c_OLS and set that to -1 or b(v))
            h_vector = zeros(1, size(c_OLS_coeff, 1));
            % Left side
            LH = LH_part(j, :);
            LH(single_nonzero_loc_1) = order; % Set to order at single-nonzero-location
            LH_idx = find(ismember(c_OLS_coeff, [1, LH], 'rows'));
            h_vector(LH_idx) = -1;
            % Right side
            RH_idx = find(ismember(c_OLS_coeff, [2*ones(size(RH, 1), 1), RH], 'rows'));
            b_v = prod(vertex_bary_21.^gamma, 2);
            h_vector(RH_idx) = b_v;
            
            % Fill H
            H = vertcat(H, h_vector);
            
        end        
    end
end

%% Equality constrained OLS estimation
% TODO: Use efficient iterative solver
Lagrangian = pinv([global_B' * global_B, H';...
            H, zeros(size(H, 1), size(H, 1))]);
C1 = Lagrangian(1:size(H, 2), 1:size(H, 2));
c_OLS = C1 * global_B' * Y_id;

% % Iterative solver
% epsilon = 10^(-6);
% lambda = ones(size(H, 1), 1);
% c_OLS = inv(2 * global_B' * global_B + 1/epsilon * H' * H) * (2 * global_B' * Y_id - H' * lambda);
% dif = 1;
% for a = 1:10
%     
%     c_OLS_new = inv(2 * global_B' * global_B + 1/epsilon * H' * H) * 2 * global_B' * global_B * c_OLS;
%     dif = abs(mean(mean(c_OLS_new - c_OLS)));
%     c_OLS = c_OLS_new;
%     
% end

global_B_val = [];
for i = 1:size(tri.ConnectivityList, 1)
    
    % Get data points and their barycentric coords in current triangle
    indices = find(IMap_val == i);
    BaryC_current = BaryC_val(indices, :);
    
    % Create sorted B-form regression matrix and global B-form matrix
    exponentials    = gen_exp(3, spline_poly_order);
    sorted_B_val    = x2fx(BaryC_current, exponentials);
    global_B_val    = blkdiag(global_B_val, sorted_B_val);
          
end


Y_est = global_B_val * c_OLS;
Y_est(Y_est > 0) = 0;
Y_est(Y_est < -0.2) = -0.2;
% % residuals = Y_est - Y_val;

%% Plotting
OLS_plotting(X_val, Y_val, Y_est, 0)

end