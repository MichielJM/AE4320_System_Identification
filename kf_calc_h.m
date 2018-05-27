%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zpred = kf_calcHx(x) Calculates the output dynamics equation h(x,u,t) 
%   
%   Author: M.J. Mollema (adapted from original by C.C. de Visser, Delft
%   University of Technology, 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zpred = kf_calcHx(t, x, u)
    
    zpred = zeros(3, 1);
    
    % output dynamics go here!
    zpred(1) = atan(x(1) / x(2)) * (1 + x(4));
    zpred(2) = atan(x(2) / sqrt(x(1)^2 + x(2)^2));
    zpred(3) = sqrt(x(1)^2 + x(2)^2 + x(3)^2);