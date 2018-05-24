%F16_IOmapping.m
%
%   Description:
%       This file shows how load and plot the IO mapping of the F16 aerodynamic coefficient C_m.
%
%   Construction date:  17-06-2009
%   Last updated:       23-03-2015
%
%   C.C. de Visser
%   TUDelft, Faculty of Aerospace Engineering, Division of Control &
%   Simulation
%
%   E. de Weerdt
%   TUDelft, Faculty of Aerospace Engineering, ASTI & Control and
%   Simulation Division
%

%%   Clearing workspace
clear all
close all
clc

printfigs = 1; % 1: print figures to disk as PNG
figpath = ''; % set the path where the figures will be printed

warning('WARNING: This file was written for Matlab R2017A, some functions may not be compatible in later versions');

%%   Loading data
load_f16data2018

%%   Plotting results
%---------------------------------------------------------
% creating triangulation (only used for plotting here)
TRIeval = delaunayn(Z_k(:, [1 2]));

%   viewing angles
az = 140;
el = 36;

% create figures

plotID = 101;
figure(plotID);
set(plotID, 'Position', [0 100 900 500], 'defaultaxesfontsize', 10, 'defaulttextfontsize', 10, 'color', [1 1 1], 'PaperPositionMode', 'auto');
% plot data points
plot3(alpha_m, beta_m, Cm, '.k'); % note that alpha_m = alpha, beta_m = beta, y = Cm
view(0, 90); 
ylabel('beta [rad]');
xlabel('alpha [rad]');
zlabel('C_m [-]');
title('F16 CM(\alpha_m, \beta_m) raw datapoints only');
% print results to disk if printfigs = 1
if (printfigs == 1)
    fpath = sprintf('fig_F16data3D');
    savefname = strcat(figpath, fpath);
    print(plotID, '-dpng', '-r300', savefname);
    fprintf('Printed figure to <%s>\n', savefname);
end



plotID = 1001;
figure(plotID);
set(plotID, 'Position', [800 100 900 500], 'defaultaxesfontsize', 10, 'defaulttextfontsize', 10, 'color', [1 1 1], 'PaperPositionMode', 'auto');
trisurf(TRIeval, alpha_m, beta_m, Cm, 'EdgeColor', 'none'); 
grid on;
hold on;
% plot data points
plot3(alpha_m, beta_m, Cm, '.k'); % note that alpha_m = alpha, beta_m = beta, y = Cm
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



