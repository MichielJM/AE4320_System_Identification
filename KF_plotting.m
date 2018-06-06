%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot relevant data from KF
%
% Author: M.J. Mollema
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = KF_plotting(printfigs, Z_k, z_pred, z_pred_corr, t)
%% Settings
figpath = 'data/'; % set the path where the figures will be printed

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

% Alpha vs Beta
plotID = 101;
figure(plotID);
set(plotID, 'Position', [150 150 720 800], 'defaultaxesfontsize', 10, 'defaulttextfontsize', 10, 'color', [1 1 1], 'PaperPositionMode', 'auto');
plot(alpha_m, beta_m, 'r', alpha_KF, beta_KF, 'b', alpha_corr, beta_KF, 'g');
view(0, 90); 
ylabel('beta [rad]');
xlabel('alpha [rad]');
title('F16 CM(\alpha_m, \beta_m)');
% print results to disk if printfigs = 1
if (printfigs == 1)
    fpath = sprintf('fig_alpha_vs_beta');
    savefname = strcat(figpath, fpath);
    print(plotID, '-dpng', '-r300', savefname);
    fprintf('Printed figure to <%s>\n', savefname);
end


% Alpha vs time
subplot(3, 1, 1); hold on;
plot(t, alpha_m, 'r')
plot(t, alpha_KF, 'b')
plot(t, alpha_corr, 'k')
xlabel('Time [sec]');
ylabel('Alpha [rad]');

% Beta vs time
subplot(3, 1, 2); hold on;
plot(t, beta_m, 'r');
plot(t, beta_KF, 'b');
xlabel('Time [sec]');
ylabel('Beta [rad]');

% V_tot vs time
subplot(3 ,1, 3); hold on;
plot(t, Vtot_m, 'r');
plot(t, Vtot_KF, 'b');
xlabel('Time [sec]');
ylabel('Vtot [m/s]');

% plotID = 1001;
% figure(plotID);
% set(plotID, 'Position', [800 100 900 500], 'defaultaxesfontsize', 10, 'defaulttextfontsize', 10, 'color', [1 1 1], 'PaperPositionMode', 'auto');
% trisurf(TRIeval, alpha_m, beta_m, Cm, 'EdgeColor', 'none'); 
% grid on;
% hold on;
% % plot data points
% plot3(alpha_m, beta_m, Cm, '.k'); % note that alpha_m = alpha, beta_m = beta, y = Cm
% view(az, el);
% ylabel('beta [rad]');
% xlabel('alpha [rad]');
% zlabel('C_m [-]');
% title('F16 CM(\alpha_m, \beta_m) raw interpolation');
% % set fancy options for plotting 
% set(gcf,'Renderer','OpenGL');
% hold on;
% poslight = light('Position',[0.5 .5 15],'Style','local');
% hlight = camlight('headlight');
% material([.3 .8 .9 25]);
% minz = min(Cm);
% shading interp;
% lighting phong;
% drawnow();
% % print results to disk if printfigs = 1
% if (printfigs == 1)
%     fpath = sprintf('fig_F16data3DSurf');
%     savefname = strcat(figpath, fpath);
%     print(plotID, '-dpng', '-r300', savefname);
%     fprintf('Printed figure to <%s>\n', savefname);
% end


