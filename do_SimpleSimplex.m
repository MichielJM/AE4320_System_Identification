% TODO: model validation, identify (at least) 5 models, clean up code

function do_SimpleSimplex(X, Y, order)

% Split data into identification and validation
X_id    = X(2:2:end, 1:2); % Only alpha and beta (???)
X_val   = X(1:2:end, 1:2);
Y_id    = Y(2:2:end);
Y_val   = Y(1:2:end);

order = 20;

% Define vertices
V_x = [1.5; -0.2; -0.2];
V_y = [0; 0.5; -0.5];
V   = [V_x, V_y];

% Find points inside simplex and get Barycentric coordinates
TRI             = delaunayTriangulation(V);
[IMap, BaryC]   = tsearchn(V, TRI, X_id);

% figure;
% plot(V_x, V_y, X(:, 1), X(:, 2))

% Created sorted B-form regression matrix
exponentials    = gen_exp(3, order);
sorted_B        = x2fx(BaryC, exponentials);

% Perform OLS to determine B-coefficients
c_OLS = pinv(sorted_B' * sorted_B) * sorted_B' * Y_id;

%% Model evaluation
[IMap, BaryC_val]   = tsearchn(V, TRI, X_val);
sorted_B_val        = x2fx(BaryC_val, exponentials);

Y_est = sorted_B_val * c_OLS;
residuals = Y_est - Y_val;


%% Plotting
% OLS_plotting(X_val, Y_val, Y_est, 1)
end