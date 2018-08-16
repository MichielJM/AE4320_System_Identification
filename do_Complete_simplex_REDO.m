function do_Complete_simplex_REDO(X, Y, spline_poly_order, spline_cont_order)

order = 0;

% Split data into identification and validation
X_id    = X(2:2:end, 1:2); % Only alpha and beta (???)
X_val   = X(1:2:end, 1:2);
Y_id    = Y(2:2:end);
Y_val   = Y(1:2:end);
% X_id = [0.2, 0.8;...
%         0.4, 0.2;...
%         0.8, 0.2;...
%         0.8, 0.6];
% Y_id = [-1; 2; 2; 1];
% X_val = X_id;
% Y_val = Y_id;

%% Define triangular grid and plot b-coefficient and vertex locations
grid_start_x    = 0; %-0.2;
grid_end_x      = 1; %0.8;
grid_start_y    = 0; %-0.3;
grid_end_y      = 1; %0.3;
% grid_start_x    = -0.2;
% grid_end_x      = 0.8;
% grid_start_y    = -0.3;
% grid_end_y      = 0.3;
step_x          = (grid_end_x - grid_start_x)/2^order;
step_y          = (grid_end_y - grid_start_y)/2^order;
[x, y]  = meshgrid(grid_start_x : step_x : grid_end_x,...
                grid_start_y : step_y : grid_end_y);
tri     = delaunayTriangulation(x(:), y(:));
trimesh(tri, x, y);
% hold on
% plot(X_id(:,1), X_id(:,2), 'x')

% Add labels to vertices in plot
vertices = tri.Points;
for i = 1:size(vertices, 1)
    vertex_label = (['v_', num2str(i - 0)]);
    text(vertices(i, 1), vertices(i, 2), vertex_label, 'Color', 'black', 'FontSize', 14);
end

% Add labels to triangles in plot
triangles = sort(tri.ConnectivityList, 2);
for i = 1:size(triangles, 1)
    triangle_label = (['t_', num2str(i)]);
    triangle_centroid = [mean(tri.Points(triangles(i, :), 1)), mean(tri.Points(triangles(i, :), 2))];
    text(triangle_centroid(1), triangle_centroid(2), triangle_label, 'Color', 'black', 'FontSize', 14)
end

% Plot B-coefficients and labels
figure
trimesh(tri, x, y)
hold on
b_cart = [];
for i = 1:size(triangles, 1)
    
    multi_index     = gen_exp(3, spline_poly_order);
    b_bary          = multi_index / spline_poly_order;
    simplex_coords  = tri.Points(triangles(i, :), :);
    b_cart          = vertcat([b_cart; bsplinen_bary2cart(simplex_coords, b_bary)]);
    
    for j = 1:size(multi_index, 1)
        b_label{j + (i-1)*size(multi_index, 1)} = (['c_{', num2str(multi_index(j, :)), '}^{t_', num2str(i), '}']);
    end

end
plot(b_cart(:, 1), b_cart(:, 2), '.g', 'MarkerSize', 25)
textfit(b_cart(:, 1), b_cart(:, 2), b_label)

%% Create global B matrix
% Find which data points are in which triangle
[IMap, BaryC] = tsearchn([x(:), y(:)], tri, X_id);
[IMap_val, BaryC_val] = tsearchn([x(:), y(:)], tri, X_val);

% Loop over all triangles to determine B-form regression matrix
c_OLS_coeff = [];
global_B = [];
global_indices = [];
for i = 1:size(triangles, 1)
    
    % Get data points and their barycentric coords in current triangle
    indices = find(IMap == i);
    global_indices = vertcat(global_indices, indices);
    BaryC_current = BaryC(indices, :);
    
    % Create sorted B-form regression matrix and global B-form matrix
    exponentials    = gen_exp(3, spline_poly_order);
    coefficients    = factorial(spline_poly_order) ./ prod(factorial(exponentials), 2);
    sorted_B        = x2fx(BaryC_current, exponentials) .* coefficients';
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
H = [];

for i = 1:size(int_edges, 1)
    
    for order = 0:spline_cont_order
        
        % Determine adjacent triangles of current edge and Out-Of-Edge
        % vertices.
        edge = int_edges(i, :);
        current_triangles = find(sum((triangles == edge(1)) + (triangles == edge(2)), 2) == 2);
        OOE_vertex_1 = triangles(current_triangles(1), :);
        OOE_vertex_1((OOE_vertex_1 == edge(1)) | (OOE_vertex_1 == edge(2))) = [];
        OOE_vertex_2 = triangles(current_triangles(2), :);
        OOE_vertex_2((OOE_vertex_2 == edge(1)) | (OOE_vertex_2 == edge(2))) = [];
        
        % Find barycentric coords of OOE_1 wrt triangle 2
        simplex_2 = tri.Points(triangles(current_triangles(2), :), :)
        OOE_bary_12 = bsplinen_cart2bary(simplex_2, tri.Points(OOE_vertex_1, :))
        
        % Find location of single non-zero value in vertex multi-index
        single_nonzero_loc_1 = find(triangles(current_triangles(1), :) == OOE_vertex_1)
        single_nonzero_loc_2 = find(triangles(current_triangles(2), :) == OOE_vertex_2)
        
        multi_index = gen_exp(3, spline_poly_order);
        indices     = exponentials(:, single_nonzero_loc_1) == order;
        LH_part     = exponentials(indices, :)
        
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
            RH_part = vertcat(RH_part, RH)
            
            % Smoothness conditions (find index of corresponding entry in
            % c_OLS and set that to -1 or b(v))
            h_vector = zeros(1, size(c_OLS_coeff, 1));
            % Left side
            LH = LH_part(j, :);
            LH(single_nonzero_loc_1) = order; % Set to order at single-nonzero-location
            LH_idx = find(ismember(c_OLS_coeff, [current_triangles(1), LH], 'rows'));
            h_vector(LH_idx) = -1;
            % Right side
            RH_idx = find(ismember(c_OLS_coeff, [current_triangles(2)*ones(size(RH, 1), 1), RH], 'rows'));
            b_v = prod(OOE_bary_12.^gamma, 2);
            h_vector(RH_idx) = b_v;
            
            % Fill H
            H = vertcat(H, h_vector);
            
        end
            
    end
    
end

end