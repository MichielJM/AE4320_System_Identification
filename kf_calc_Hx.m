%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H = kf_calcDHx(x) Calculates the Jacobian of the output dynamics equation h(x,u,t) 
%   
%   Author: M.J. Mollema (adapted from original by C.C. de Visser, Delft
%   University of Technology, 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hx = kf_calcDHx(t, x, u)

    Hx = zeros(length(u), length(x));
    
    % Calculate Jacobian matrix of output dynamics
    H00 = 1/(x(2) + x(1)^2/x(2)^3) * (1 + x(4));
    H01 = -x(1)/(x(2)^2 + x(1)^2/x(2)^4) * (1 + x(4));
    H02 = 0;
    H03 = atan(x(1)/x(2));
    
    H10 = (x(1) * x(2)) / ((x(1)^2 + x(3)^2)^(3/2) + x(2)^2 / sqrt(x(1)^2 + x(3)^2));
    H11 = 1 / (sqrt(x(1)^2 + x(3)^2) + x(2)^2 / sqrt(x(1)^2 + x(3)^2));
    H12 = (x(3) * x(2)) / ((x(1)^2 + x(3)^2)^(3/2) + x(2)^2 / sqrt(x(1)^2 + x(3)^2));
    H13 = 0;
    
    H20 = x(1) / sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    H21 = x(2) / sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    H22 = x(3) / sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    H23 = 0;
    
    Hx = [H00, H01, H02, H03;...
        H10, H11, H12, H13;...
        H20, H21, H22, H23];