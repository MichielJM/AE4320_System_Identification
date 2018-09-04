function do_OLS(X, Y, interpolation_order)
% DO_OLS Performs an ordinary least squares estimation on a given dataset.
% 
% Inputs:
%  X: Data-points to estimate with
%  Y: Data-points to estimate
%  interpolation_order: order of the estimating polynomial
% 
% Output:
%  - None
% 
% M.J. Mollema - 04.09.2018

% Split data into identification and validation set (odd index for
% validation, even for identification)
Y_id   = Y(2:2:end);
Y_val  = Y(1:2:end);
X_id    = X(2:2:end, :);
X_val   = X(1:2:end, :);

% Get theta
Ax = create_regression_matrix(X_id(:, 1:3), interpolation_order);
theta_OLS = pinv(Ax' * Ax) * Ax' * Y_id;

% Get Ax for validation data
Ax_val = create_regression_matrix(X_val(:, 1:3), interpolation_order);

% Estimate Cm for validation data
Cm_est_val = Ax_val * theta_OLS;

%%   Plotting OLS
% if showfigs
%     OLS_plotting(X_val, Y_val, Cm_est_val, printfigs)
% end

%% Model validation
% Statistical (use covariances)
Cov = cov(theta_OLS);

% Model-error based
residual = Y_val - Cm_est_val;
res_RMS = sqrt(mean(residual.^2));

conf = 1.96 / sqrt(length(residual));
[acx,lags] = xcorr(residual-mean(residual), 'coeff');
k_0 = find(lags == 0);
corr_test = abs(acx(k_0 + 1 : end)) / acx(k_0);

% Determine "whiteness"
inside_bounds = length(find(corr_test <= conf));
perc_inside = inside_bounds / length(corr_test) * 100

% Plot model error autocorrelation
% figure; hold on;
% line([lags(1), lags(end)], [conf, conf], 'Color','red','LineStyle','--')
% line([lags(1), lags(end)], [-conf, -conf], 'Color','red','LineStyle','--')
% plot(lags, acx)

end