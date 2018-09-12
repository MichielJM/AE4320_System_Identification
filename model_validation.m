function [residual, acx, lags, VAR, conf] = model_validation(Y_val, Y_est, Ax)
% MODEL_VALIDATION Performs residual based and statistical model
% validation.
% 
% Inputs:
%  - Y_val: validation Y-values
%  - Y_est: estimated Y-values
%  - Ax: regression matrix
% 
% Outputs:
%  - residual: model residual
%  - acx: residual autocorrelation
%  - lags: lags for autocorrelation
%  - VAR: B-coefficient variances
% (max(Y_val) - min(Y_val))
% 
% M.J. Mollema - 09.09.2018

%% Model-residual based
% For OLS estimator to be BLUE: E(residuals) = 0, residuals uncorrelated
residual = Y_val - Y_est;
res_RMS = rms(residual);
res_RMS_rel = res_RMS / (max(Y_val) - min(Y_val));
res_rel_max = max(residual / (max(Y_val) - min(Y_val)));

% Check correlation
conf = 1.96 / sqrt(length(residual));
[acx,lags] = xcorr(residual-mean(residual), 'coeff');
perc_inside_bounds = length(acx(abs(acx) < conf)) / length(acx) * 100;

% Display results
fprintf("Relative RMS of residuals: %.3f %%\n", res_RMS_rel * 100)
fprintf("Maximum relative residual: %.2f %%\n", res_rel_max*100)
fprintf("Mean of residual: %.3f \n", mean(residual))
fprintf("Percentage inside 95%% confidence bounds: %.1f \n", perc_inside_bounds)

%% Statistical validation (variances of coefficients)
COV = pinv(Ax' * Ax);
sigma = (residual' * residual) / (size(Ax, 1) - size(Ax, 2));
VAR = sigma * diag(COV);

end