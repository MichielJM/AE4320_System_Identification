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

% Define triangular grid
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
    vertex_label = (['v_', num2str(i - 1)]);
    text(vertices(i, 1), vertices(i, 2), vertex_label, 'Color', 'black', 'FontSize', 14);
end

% Add labels to triangles in plot
triangles = sort(tri.ConnectivityList, 2);
for i = 1:size(triangles, 1)
    triangle_label = (['t_', num2str(i)]);
    triangle_centroid = [mean(tri.Points(triangles(i, :), 1)), mean(tri.Points(triangles(i, :), 2))];
    text(triangle_centroid(1), triangle_centroid(2), triangle_label, 'Color', 'black', 'FontSize', 14)
end

% Plot B-coefficients and labels REFERENCE VERTEX IS FIRST IN ARRAY


figure
trimesh(tri, x, y)
hold on
for i = 1:size(triangles, 1)
    multi_index     = gen_exp(3, spline_poly_order);
    b_bary          = multi_index / spline_poly_order;
    simplex_coords  = tri.Points(triangles(i, :), :);
    b_cart          = bsplinen_bary2cart(simplex_coords, b_bary);
    
    for j = 1:size(multi_index, 1)
        b_label{j} = (['c_{', num2str(multi_index(j, :)), '}^{t_', num2str(i), '}']);
    end

    
    if mod(i, 2) == 0
        VA = 'top';
    else
        VA = 'bottom';
    end
    plot(b_cart(:, 1), b_cart(:, 2), '.g', 'MarkerSize', 25)
    text(b_cart(:, 1), b_cart(:, 2), b_label, 'VerticalAlignment', VA)
end


end