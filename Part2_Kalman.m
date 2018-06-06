%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2 of the AE4320:System Identification assignment
% A Kalman filter is implemented to estimate a bias in AoA
% measurements in the given dataset
%
%   Author: M.J. Mollema (adapted from original by C.C. de Visser, Delft)
%   University of Technology, 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all, clc, clear all

%% Load data
filename = 'data/F16traindata_CMabV_2018';
load(filename, 'Cm', 'Z_k', 'U_k');

% transpose
Cm = Cm'; Z_k = Z_k'; U_k = U_k';

%% Set simulation parameters
dt              = 0.01;
t               = 0:dt:dt*(size(Z_k, 2)-1);

stdv    = [0.01, 0.0058, 0.112];
stdw    = [1e-3, 1e-3, 1e-3, 0];

%% Run IEKF
printfigs = 0;

% Check observability of the states
check_observability

[XX_k1k1, z_pred, IEKFitcount] = IEKF(U_k, Z_k, dt, stdv, stdw);

% Correct alpha for bias using estimate bias
z_pred_corr = z_pred;
z_pred_corr(1, :) = z_pred_corr(1, :) ./ (1 + XX_k1k1(4, :));

%% Plotting
F16_PlotData
