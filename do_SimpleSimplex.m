function do_SimpleSimplex(X, Y, order, plot)
% DO_SIMPLESIMPLEX Performs parameter estimation using the simplex method
% on a single simplex.
% 
% Inputs:
%  - X: data to estimate with
%  - Y: data to estimate
%  - order: order of the estimating polynomial
%  - plot: Boolean to determine whether to plot or not
% 
% Output:
%  - None
% 
% M.J.Mollema - 07.09.2018

do_coeff = 1;
% Split data into identification and validation
[X_id, Y_id, X_val, Y_val] = split_data(X, Y, 0.5);
% X_id = X_id'; Y_id = Y_id'; X_val = X_val'; Y_val = Y_val';

% Define vertices
V_x = [1.5; -0.2; -0.2];
V_y = [0; 0.5; -0.5];
V   = [V_x, V_y];

% Find points inside simplex and get Barycentric coordinates
TRI          = delaunayTriangulation(V);
[~, BaryC]   = tsearchn(V, TRI, X_id);

% Created sorted B-form regression matrix
exponentials    = gen_exp(3, order);
coefficients    = factorial(order) ./ prod(factorial(exponentials), 2); % not needed, OLS itself sorts it out anyways
if do_coeff
    sorted_B        = x2fx(BaryC, exponentials) .* coefficients';
else
    sorted_B        = x2fx(BaryC, exponentials);
end

% Perform OLS to determine B-coefficients
c_OLS = pinv(sorted_B' * sorted_B) * sorted_B' * Y_id;

%% Model evaluation
[~, BaryC_val]   = tsearchn(V, TRI, X_val);
sorted_B_val        = x2fx(BaryC_val, exponentials);

if do_coeff
    Y_est = sorted_B_val .* coefficients' * c_OLS;
else
    Y_est = sorted_B_val * c_OLS;
end
residuals = Y_est - Y_val;


%% Plotting
if plot
    figure; hold on;
    TRI_eval = delaunayn(X_val);
    grey = [192, 192, 192]/255;
    trisurf(TRI_eval, X_val(:, 1), X_val(:, 2), Y_val,...
        'Facecolor', grey, 'FaceAlpha', 0.8, 'EdgeColor', 'None')
    plot3(X_val(:, 1), X_val(:, 2), Y_est, '.g');
end

end