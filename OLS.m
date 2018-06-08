% Function to perform Ordinary Least Squares estimation
%
% Inputs:
%   X = variables to estimate with (N_variables, N_datapoints)
%   Y = variable to estimate (1, N_datapoints)
%   order = order of estimator
%
% Outputs:
%   Ax = regression matrix (N_datapoints, order+1)
%   theta = parameter vector (order+1, 1)
%
% M.J. Mollema, 08/06/2018

function [Ax, theta] = OLS(X, Y, order)

Ax = create_regression_matrix(X(1:2, :), order);

theta = pinv(Ax' * Ax) * Ax' * Y'; %Y' to get size (datapoints, 1)
end