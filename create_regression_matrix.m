function A = create_regression_matrix(X, order)
% CREATE_REGRESSION_MATRIX Creates a regression matrix of a certain order,
%   including all cross-terms.
%
% Inputs:
% - X: state vector, shape (N_states, N)
% - order: order of the polynomial
%
% Outputs:
% - A: regression matrix
%
% M.J. Mollema (adapted from: Jesse Hagenaars - 06.05.2018)

N_states = size(X, 2);

exponentials = zeros(1, N_states);

% Go over all orders, and stack
for d = 1:order
    
    % Generate all possible combinations that sum op to order d
    exponentials = [exponentials; gen_exp(N_states, d)];
    
end

% Create regression matrix using MATLAB's x2fx
A = x2fx(X, exponentials);

end
