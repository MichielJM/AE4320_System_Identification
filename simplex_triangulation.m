function [tri, x, y] = simplex_triangulation(num_triangles_x, num_triangles_y,...
    spline_poly_order, plot_triangulation)
% SIMPLEX_TRIANGULATION Create an even-spaced triangulation for the data
% and generate the locations and names of all vertices and b-coefficients.
%
% Inputs:
%  - num_triangles_x: desired number of triangles in x-direction
%  - num_triangles_y: desired number of triangles in y-direction
%  - spline_poly_order: order of the polynomial spline, used to determine
%     the amount of and location of b-coefficients
%  - plot_triangulation: Boolean to plot or not
%     
% Outputs:
% - tri: MATLAB triangulation object
% - x: x coordinates of meshgrid of vertices
% - y: y coordinates of meshgrid of vertices
% - triangles: array containing different triangles in each row and their 
%     corresponding vertices in the columns.
% 
% M.J. Mollema - 03.09.2018

%% Define outer edges to encompass all the data with triangulation
grid_start_x    = -0.2;
grid_end_x      = 0.8;
grid_start_y    = -0.3;
grid_end_y      = 0.3;
step_x          = (grid_end_x - grid_start_x)/num_triangles_x;
step_y          = (grid_end_y - grid_start_y)/num_triangles_y;
[x, y]  = meshgrid(grid_start_x : step_x : grid_end_x,...
                grid_start_y : step_y : grid_end_y);
            
% Create triangulation
tri     = delaunayTriangulation(x(:), y(:));
triangles = sort(tri.ConnectivityList, 2);

%% Plotting
if plot_triangulation
    
    trimesh(tri, x, y);

    % Add labels to vertices in plot
    vertices = tri.Points;
    for i = 1:size(vertices, 1)
        vertex_label = (['v_', num2str(i - 0)]);
        text(vertices(i, 1), vertices(i, 2), vertex_label, 'Color', 'black', 'FontSize', 14);
    end

    % Add labels to triangles in plot
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

end
    
end