%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zpred = kf_calcHx(x) Calculates the output dynamics equation h(x,u,t) 
%   
%   Author: M.J. Mollema (adapted from original by C.C. de Visser, Delft
%   University of Technology, 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zpred = kf_calcHx(t, x, U)
    
    u = x(1); v = x(2); w = x(3); C = x(4);

    zpred = zeros(3, 1);
    
    % output dynamics go here!
    zpred(1) = atan(w / u) * (1 + C);
    zpred(2) = atan(v / sqrt(u^2 + w^2));
    zpred(3) = sqrt(u^2 + v^2 + w^2);
