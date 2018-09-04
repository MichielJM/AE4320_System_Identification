function [X, Y] = do_Kalman(Z_k, U_k, stdv, stdw, plot_kalman)
% DO_KALMAN Perform Kalman filtering on the provided set of data
% 
% Inputs:
%  - Z_k: measurement vector
%  - U_k: input vector
%  - stdv: sensor noise statistics
%  - stdw: process noise statistics
%  - plot_kalman: Boolean determining to plot or not
% 
% Outputs:
%  - X: Kalman filtered measurements
% 
% M.J. Mollema - 03.09.2018

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
% Y = Cm';

%% Plotting KF
% if plot_kalman
%     KF_plotting(printfigs, Z_k, z_pred, z_pred_corr, t)
% end