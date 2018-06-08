function exponentials = gen_exp(N_states, order)
% GEN_EXP Generate all possible combinations of exponentials that sum up to
%   a certain order.
%
% Inputs:
% - N_states: number of states
% - order: order to sum up to
%
% Outputs:
% - exponentials: matrix with all possible combinations of exponentials
%   (that sum up to order) in its rows
%
% Jesse Hagenaars - 06.05.2018

if N_states <= 1
    
    exponentials = order;
    
else
    
    % Create new row
    exponentials = zeros(0, N_states);
    
    % Go over all remaining orders
    for i = order:-1:0
        
        % Recursively call function, decreasing order by 1
        rc = gen_exp(N_states - 1, order - i);
        exponentials = [exponentials; i * ones(size(rc,1), 1), rc];

    end

end