%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H = kf_calcDHx(x) Calculates the Jacobian of the output dynamics equation h(x,u,t) 
%   
%   Author: M.J. Mollema (adapted from original by C.C. de Visser, Delft
%   University of Technology, 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hx = kf_calcDHx(t, x, U)

    u = x(1); v = x(2); w = x(3); C = x(4);
    
    Hx = zeros(length(U), length(x));
    
    % Calculate Jacobian matrix of output dynamics
    H00 = -(C + 1) * w / (u^2 + w^2);
    H01 = 0;
    H02 = (C + 1) * u / (u^2 + w^2);
    H03 = atan(w/u);
    
    H10 = -(u * v) / (sqrt(u^2 + w^2) * (u^2 + v^2 + w^2));
    H11 = sqrt(u^2 + w^2) / (u^2 + v^2 + w^2);
    H12 = -(v * w) / (sqrt(u^2 + w^2) * (u^2 + v^2 + w^2));
    H13 = 0;
    
    H20 = u / sqrt(u^2 + v^2 + w^2);
    H21 = v / sqrt(u^2 + v^2 + w^2);
    H22 = w / sqrt(u^2 + v^2 + w^2);
    H23 = 0;
    
    Hx = [H00, H01, H02, H03;...
        H10, H11, H12, H13;...
        H20, H21, H22, H23];