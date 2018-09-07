function [X_id, Y_id, X_val, Y_val] = split_data(X, Y, perc_val)
% SPLIT_DATA Splits the dataset into a validation and identification set
% with random permutations.
% 
% Inputs:
%  - X: X-values from dataset
%  - Y: Y-values from dataset
%  - perc_val: percentage of dataset to use for validation
% 
% Outputs:
%  - X_id: X-values for identification
%  - Y_id: Y-values for identification
%  - X_val: X-values for validation
%  - Y_val: Y-values for validation
% 
% M.J. Mollema - 07.09.2018

% Get random permutations
p       = randperm(length(X(:, 1)));
p_val   = p(1 : floor(perc_val*length(X(:, 1))));
p_id    = p(floor(perc_val*length(X(:, 1)))+1 : end);

% Split data
X_id    = X(p_id, 1:2);
Y_id    = Y(p_id);
X_val   = X(p_val, 1:2);
Y_val   = Y(p_val);

end