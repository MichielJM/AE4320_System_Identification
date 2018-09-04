%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F = kf_calcDFx(x) Calculates the Jacobian of the system dynamics equation f(x,u,t) 
%   
%   Author: M.J. Mollema (adapted from original by: C.C. de Visser, Delft
%   University of Technology, 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DFx = kf_calc_Fx(t, x, u)
% KF_CALC_FX Calculates the Jacobian of the system dynamics equation
% f(x,u,t) for use in the Kalman filter.
% 
% Inputs:
%  - t: time
%  - x: one-step-ahead state prediction
%  - u: input vector
% 
% Output:
%  DFx: Jacobian of system dynamics equation
% 
% M.J. Mollema (adapted from C.C. de Visser, Delft) - 04.09.2018

    % Calculate Jacobian matrix of system dynamics
    DFx = zeros(4);
    
end