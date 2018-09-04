function zpred = kf_calc_h(t, x, U)
% KF_CALC_H Calculates the output dynamics equation h(x,u,t) for use in the
% Kalman filter.
% 
% Input:
%  - t: time
%  - x: one-step-ahead optimal state estimation
%  - u: input vector
% 
% Output:
%  - zpred: output dynamics equations
% 
% M.J. Mollema (adapted from C.C. de Visser, Delft) - 04.09.2018
    
    u = x(1); v = x(2); w = x(3); C = x(4);

    zpred = zeros(3, 1);
    
    % output dynamics go here!
    zpred(1) = atan(w / u) * (1 + C);
    zpred(2) = atan(v / sqrt(u^2 + w^2));
    zpred(3) = sqrt(u^2 + v^2 + w^2);

end