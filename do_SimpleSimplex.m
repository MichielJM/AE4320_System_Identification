function do_SimpleSimplex(X, Y, order, plotting, plot_triangulation)
% DO_SIMPLESIMPLEX Performs parameter estimation using the simplex method
% on a single simplex.
% 
% Inputs:
%  - X: data to estimate with
%  - Y: data to estimate
%  - order: order of the estimating polynomial
%  - plotting: Boolean to determine whether to plot or not
% 
% Output:
%  - None
% 
% M.J.Mollema - 07.09.2018

fprintf("\nPerforming single simplex estimation with order %.f polynomial\n", order)

% Split data into identification and validation
[X_id, Y_id, X_val, Y_val] = split_data(X, Y, 0.5);

% Define vertices
V_x = [1.5; -0.2; -0.2];
V_y = [0; 0.5; -0.5];
V   = [V_x, V_y];

% Find points inside simplex and get Barycentric coordinates
tri          = delaunayTriangulation(V);
triangles    = sort(tri.ConnectivityList, 2);
[~, BaryC]   = tsearchn(V, tri, X_id);

% Created sorted B-form regression matrix
exponentials    = gen_exp(3, order);
sorted_B        = x2fx(BaryC, exponentials);

% Perform OLS to determine B-coefficients
c_OLS = pinv(sorted_B' * sorted_B) * sorted_B' * Y_id;

%% Model evaluation
[~, BaryC_val]   = tsearchn(V, tri, X_val);
sorted_B_val     = x2fx(BaryC_val, exponentials);
Y_est = sorted_B_val * c_OLS;

[residual, acx, lags, VAR, conf] = model_validation(Y_val, Y_est, sorted_B);

%% Plotting
if plotting
    % Plot estimation result
    figure; hold on;
    view(20, 12);
    TRI_eval = delaunayn(X_val);
    grey = [192, 192, 192]/255;
    trisurf(TRI_eval, X_val(:, 1), X_val(:, 2), Y_val,...
        'Facecolor', grey, 'FaceAlpha', 0.8, 'EdgeColor', 'None')
    plot3(X_val(:, 1), X_val(:, 2), Y_est, '.g');
    
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

if plot_triangulation
    
    font_size = 16;
    figure; hold on;
    trimesh(tri, V_x, V_y);
    plot(X(:, 1), X(:, 2));
    xlabel('\alpha [rad]', 'fontsize', font_size)
    ylabel('\beta [rad]', 'fontsize', font_size)

    % Add labels to vertices in plot
    vertices = tri.Points;
    for i = 1:size(vertices, 1)
        vertex_label = (['v_', num2str(i - 0)]);
        text(vertices(i, 1), vertices(i, 2), vertex_label, 'Color', 'black', 'FontSize', font_size);
    end

    % Add labels to triangles in plot
    for i = 1:size(triangles, 1)
        triangle_label = (['t_', num2str(i)]);
        triangle_centroid = [mean(tri.Points(triangles(i, :), 1)), mean(tri.Points(triangles(i, :), 2))];
        text(triangle_centroid(1), triangle_centroid(2), triangle_label, 'Color', 'black', 'FontSize', font_size)
    end

    % Plot B-coefficients and labels
    figure; hold on;
    trimesh(tri, V_x, V_y)
    xlabel('\alpha [rad]', 'fontsize', font_size)
    ylabel('\beta [rad]', 'fontsize', font_size)
    b_cart = [];
    for i = 1:size(triangles, 1)

        multi_index     = gen_exp(3, order);
        b_bary          = multi_index / order;
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