% Main file for the AE4320 Assignment: Multivariate Simplex Splines.
% Running this file will go trough all the different parts of the
% assignment, from Kalman filtering to simplex splines system
% identification. The settings section below is used to set different
% interpolation order and that stuff, as well of the option to turn on/off
% different plotting parts.
% 
% M.J. Mollema - 04.09.2018

close all, clc, clear all

%% Settings
% interpolation orders
interpolation_order  = 8;
simple_simplex_order = 2;
spline_poly_order    = 10;
spline_cont_order    = 2;

% Amount of simplices
num_simplices_x = 2;
num_simplices_y = 2;

% Plotting figures
plot_kalman         = 0;
plot_OLS            = 0;
plot_simple_simplex = 0;
plot_triangulation  = 0; % When having high order or lots of simplices this plotting can take some time
plot_compl_simplex  = 0;

%% Load data and set noise statistics
filename = 'data/F16traindata_CMabV_2018';
load(filename, 'Cm', 'Z_k', 'U_k');

% transpose
Cm = Cm'; Z_k = Z_k'; U_k = U_k'; 

% Noise STD
stdv    = [0.01, 0.0058, 0.112];
stdw    = [1e-3, 1e-3, 1e-3, 0];

%% Perform Kalman filtering
[X_kalman] = do_Kalman(Z_k, U_k, stdv, stdw, plot_kalman);
Y_kalman = Cm';

%% Ordinary least squares estimator
do_OLS(X_kalman, Y_kalman, interpolation_order)

%% Perform single simplex polynomial
do_SimpleSimplex(X_kalman, Y_kalman, simple_simplex_order)

%% Perform complete simplex splines
do_CompleteSimplex(X_kalman, Y_kalman, spline_poly_order, spline_cont_order,...
    num_simplices_x, num_simplices_y, plot_triangulation)











