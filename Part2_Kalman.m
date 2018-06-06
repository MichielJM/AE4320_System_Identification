%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2 of the AE4320:System Identification assignment
% A Kalman filter is implemented to estimate a bias in AoA
% measurements in the given dataset
%
%   Author: M.J. Mollema (adapted from original by C.C. de Visser, Delft)
%   University of Technology, 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all, clc, clear all
printfigs = 0;
showfigs = 0;

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

%% Plotting
if showfigs
    KF_plotting(printfigs, Z_k, z_pred, z_pred_corr, t)
end
%% Ordinary least squares estimator
% Steps:
% 1 Obtain measurement data
% 2 Formulate linear regression model structure
% 3 Formulate regression matrix
% 4 Formulate least squares estimator
% 5 Evaluate/validate model

% Step 1
alpha = z_pred_corr(1, :);
beta = z_pred_corr(2, :);
Cm = Cm;

% Step 2: p(x,theta) = A(x) * theta
order = 1;
% order_x = 1;
% order_y = 1;

% Step 3: A(x)
sz_A = 0.5*order^2 + 1.5*order + 1;
Ax = zeros(size(z_pred, 2), sz_A);
for i = 1:size(z_pred, 2)
    Ax(i, :) = [1, alpha(i), beta(i)];
end

% Step 4: theta_OLS
theta_OLS = (Ax' * Ax)^(-1) * Ax' * Cm'; %Cm' to get size (datapoints, 1)

















