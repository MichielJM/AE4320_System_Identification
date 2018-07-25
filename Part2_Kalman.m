% Part 2 of the AE4320:System Identification assignment
% A Kalman filter is implemented to estimate a bias in AoA
% measurements in the given dataset
%
%   Author: M.J. Mollema (adapted from original by C.C. de Visser, Delft)
%   University of Technology, 2013)

close all, clc, clear all
printfigs = 0;
showfigs = 0;
interpolation_order = 3;
simple_spline_order = 2;

%% Load data and set noise statistics
filename = 'data/F16traindata_CMabV_2018';
load(filename, 'Cm', 'Z_k', 'U_k');

% transpose
Cm = Cm'; Z_k = Z_k'; U_k = U_k'; 

% Noise STD
stdv    = [0.01, 0.0058, 0.112];
stdw    = [1e-3, 1e-3, 1e-3, 0];

%% Set time vector
dt              = 0.01;
t               = 0:dt:dt*(size(Z_k, 2)-1);

%% Run IEKF
% Check observability of the states
check_observability

% Run IEKF script
[XX_k1k1, z_pred, IEKFitcount] = IEKF(U_k, Z_k, dt, stdv, stdw);

% Correct alpha for bias using estimate bias
z_pred_corr = z_pred;
z_pred_corr(1, :) = z_pred_corr(1, :) ./ (1 + XX_k1k1(4, :));

X = z_pred_corr';
Y = Cm';

%% Plotting KF
% if showfigs
%     KF_plotting(printfigs, Z_k, z_pred, z_pred_corr, t)
% end

%% Ordinary least squares estimator CHECK REGRESSION MATRIX
% Split data into identification and validation set (odd index for
% validation, even for identification)
Y_id   = Cm(2:2:end)';
Y_val  = Cm(1:2:end)';
X_id    = z_pred_corr(:, 2:2:end)';
X_val   = z_pred_corr(:, 1:2:end)';

% Shuffle
% perm = randperm(size(Cm_id, 1));
% Cm_id   = Cm_id(perm);
% z_id    = z_id(perm, :);
% 
% perm = randperm(size(Cm_val, 1));
% Cm_val  = Cm_val(perm);
% z_val   = z_val(perm, :);

% Get theta
[Ax, theta_OLS] = OLS(X_id, Y_id, interpolation_order);

% Get Ax for validation data
Ax_val = create_regression_matrix(X_val, interpolation_order);

% Estimate Cm for validation data
Cm_est_val = Ax_val * theta_OLS;

%%   Plotting OLS
if showfigs
    OLS_plotting(X_val, Y_val, Cm_est_val, printfigs)
end

%% Model validation
% Statistical (use covariances)
Cov = cov(theta_OLS)

% Model-error based
residual = Cm_est_val - Y_val;

conf = 1.96 / sqrt(length(residual));
[acx,lags] = xcorr(residual-mean(residual), 'coeff');
acx_norm = acx ./ max(acx, 1);

% Plot model error autocorrelation
% figure; hold on;
% line([lags(1), lags(end)], [conf, conf], 'Color','red','LineStyle','--')
% line([lags(1), lags(end)], [-conf, -conf], 'Color','red','LineStyle','--')
% plot(lags, acx_norm)

% Determine "whiteness"
% inside_bounds   = length(find(-conf < acx < conf));
% perc_inside     = inside_bounds / length(acx) * 100

%% Perform single simplex polynomial
do_SimpleSimplex(X, Y, simple_spline_order)











