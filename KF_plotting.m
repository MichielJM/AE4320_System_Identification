function [] = KF_plotting(Z_k, z_pred, z_pred_corr, t, IEKFitcount, C_a_up)
% KF_PLOTTING Plots the relevant things from the Kalman filter
% 
% Inputs:
%  - Z_k: original measurement vector
%  - z_pred: Kalman filtered measurement vector
%  - z_pred_corr: Kalman filtered and bias corrected measurement vector
%  - t: time vector
%  - IEKFitcount: number of iterations vector for IEKF
% 
% Output:
%  - None
% M.J. Mollema - 07.09.2018

%% Set variable names
% Measured
alpha_m = Z_k(1, :);
beta_m = Z_k(2, :);
Vtot_m = Z_k(3, :);

% Kalman filtered
alpha_KF = z_pred(1, :);
beta_KF = z_pred(2, :);
Vtot_KF = z_pred(3, :);

% Bias corrected (only alpha has bias)
alpha_corr = z_pred_corr(1, :);
beta_corr = z_pred_corr(2, :); % = beta_KF
Vtot_corr = z_pred_corr(3, :); % = Vtot_KF

%%   Plotting results
fontsize = 16;
last_datapoint = 1000;
% Alpha vs Beta
plotID = 101;
figure(plotID);
set(plotID, 'defaultaxesfontsize', fontsize, 'defaulttextfontsize', fontsize, 'color', [1 1 1], 'PaperPositionMode', 'auto');
plot(alpha_m, beta_m, 'r', alpha_KF, beta_KF, 'b', alpha_corr, beta_KF, 'g');
pbaspect([1.6/0.6 1 1])
view(0, 90); 
ylabel('\beta [rad]');
xlabel('\alpha [rad]');

% Alpha vs time
plotID = 201;
figure(plotID)
set(plotID, 'defaultaxesfontsize', fontsize, 'defaulttextfontsize', fontsize, 'color', [1 1 1], 'PaperPositionMode', 'auto');
subplot(3, 1, 1); hold on;
plot(t(1:last_datapoint), alpha_m(1:last_datapoint), 'r')
plot(t(1:last_datapoint), alpha_KF(1:last_datapoint), 'b')
plot(t(1:last_datapoint), alpha_corr(1:last_datapoint), 'g')
xlabel('Time [sec]');
ylabel('\alpha [rad]');

% Beta vs time
subplot(3, 1, 2); hold on;
plot(t(1:last_datapoint), beta_m(1:last_datapoint), 'r');
plot(t(1:last_datapoint), beta_KF(1:last_datapoint), 'b');
xlabel('Time [sec]');
ylabel('\beta [rad]');

% V_tot vs time
subplot(3 ,1, 3); hold on;
plot(t(1:last_datapoint), Vtot_m(1:last_datapoint), 'r');
plot(t(1:last_datapoint), Vtot_KF(1:last_datapoint), 'b');
xlabel('Time [sec]');
ylabel('V_{tot} [m/s]');

% IEKF iteration count
plotID = 301;
figure(plotID)
set(plotID, 'defaultaxesfontsize', fontsize, 'defaulttextfontsize', fontsize, 'color', [1 1 1], 'PaperPositionMode', 'auto');
plot(t, IEKFitcount)
xlabel('Time [sec]')
ylabel('Number of iterations [-]')

% Upwash coefficient
plotID = 401;
figure(plotID)
set(plotID, 'defaultaxesfontsize', fontsize, 'defaulttextfontsize', fontsize, 'color', [1 1 1], 'PaperPositionMode', 'auto');
plot(t, C_a_up)
xlabel('Time [sec]')
ylabel('Upwash coefficient C_{\alpha_{up}} [-]')
