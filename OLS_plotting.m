% Plotting of OLS estimator results
%
% Input:
%     Z = measurement array (N_variables, N_datapoints)
%     Y = values to estimate (1, N_datapoints)
%     Y_est = estimate values of Y (N_datapoints, 1)
%     printfigs = flag to save figs to disk or not
% Output = [] (just shows plots)
% 
% M.J. Mollema, 08/06/2018
function [] = OLS_plotting(Z, Y, Y_est, printfigs)

% Get data
alpha = Z(1, :);
beta = Z(2,:);
Vtot = Z(3, :);
Cm = Y;
Cm_est = Y_est;

% creating triangulation (only used for plotting here)
TRIeval = delaunayn(Z([1 2], :)');

%   viewing angles
az = 140;
el = 36;

plotID = 1001;
figure(plotID);
set(plotID, 'Position', [800 100 900 500], 'defaultaxesfontsize', 10, 'defaulttextfontsize', 10, 'color', [1 1 1], 'PaperPositionMode', 'auto');
trisurf(TRIeval, alpha', beta', Cm_est, 'EdgeColor', 'none'); 
grid on;
hold on;
% plot data points
plot3(alpha', beta', Cm, '.k'); % note that alpha_m = alpha, beta_m = beta, y = Cm
view(az, el);
ylabel('beta [rad]');
xlabel('alpha [rad]');
zlabel('C_m [-]');
title('F16 CM(\alpha_m, \beta_m) raw interpolation');
% set fancy options for plotting 
set(gcf,'Renderer','OpenGL');
hold on;
poslight = light('Position',[0.5 .5 15],'Style','local');
hlight = camlight('headlight');
material([.3 .8 .9 25]);
minz = min(Cm);
shading interp;
lighting phong;
drawnow();
% print results to disk if printfigs = 1
if (printfigs == 1)
    fpath = sprintf('fig_F16data3DSurf');
    savefname = strcat(figpath, fpath);
    print(plotID, '-dpng', '-r300', savefname);
    fprintf('Printed figure to <%s>\n', savefname);
end
end