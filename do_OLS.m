function do_OLS(X, Y, interpolation_order, plotting)
% DO_OLS Performs an ordinary least squares estimation on a given dataset.
% 
% Inputs:
%  - X: Data-points to estimate with
%  - Y: Data-points to estimate
%  - interpolation_order: order of the estimating polynomial
%  - plotting: Boolean to determine whether to plot or not
% 
% Output:
%  - None
% 
% M.J. Mollema - 04.09.2018

fprintf("\nPerforming OLS estimation with order %.f polynomial\n", interpolation_order)

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
[residual, acx, lags, VAR, conf] = model_validation(Y_val, Y_est, Ax);

%% Plotting
if plotting
    % Plot estimation result
    figure; hold on;
    view(20,12)
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

end