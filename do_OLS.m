function do_OLS(X, Y, interpolation_order, plot)
% DO_OLS Performs an ordinary least squares estimation on a given dataset.
% 
% Inputs:
%  - X: Data-points to estimate with
%  - Y: Data-points to estimate
%  - interpolation_order: order of the estimating polynomial
%  - plot: Boolean to determine whether to plot or not
% 
% Output:
%  - None
% 
% M.J. Mollema - 04.09.2018

% Split data into identification and validation set
[X_id, Y_id, X_val, Y_val] = split_data(X, Y, 0.5);

% Get theta
Ax = create_regression_matrix(X_id, interpolation_order);
theta_OLS = pinv(Ax' * Ax) * Ax' * Y_id;

% Get Ax for validation data
Ax_val = create_regression_matrix(X_val, interpolation_order);

% Estimate Cm for validation data
Y_est = Ax_val * theta_OLS;

%% Model validation
% Model-residual based
residual = Y_val - Y_est;
res_RMS = sqrt(mean(residual.^2));
res_RMS_rel = res_RMS / (max(Y_val) - min(Y_val));

% conf = 1.96 / sqrt(length(residual));
% [acx,lags] = xcorr(residual-mean(residual), 'coeff');
% k_0 = find(lags == 0);
% corr_test = abs(acx(k_0 + 1 : end)) / acx(k_0);
% 
% % Determine "whiteness"
% inside_bounds = length(find(corr_test <= conf));
% perc_inside = inside_bounds / length(corr_test) * 100
% 
% % Plot model error autocorrelation
% % figure; hold on;
% % line([lags(1), lags(end)], [conf, conf], 'Color','red','LineStyle','--')
% % line([lags(1), lags(end)], [-conf, -conf], 'Color','red','LineStyle','--')
% % plot(lags, acx)

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