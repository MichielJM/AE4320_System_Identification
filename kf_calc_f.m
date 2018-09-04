function xdot = kf_calc_f(t, x, u)
% KF_CALC_F Calculates the system dynamics equation f(x,u,t) for use in the
% Kalman filter.
% 
% Inputs:
%  - t: time (for integration using rk4)
%  - x: one-step-ahead optimal state estimation
%  - u: input vector
% 
% Output:
%  xdot: system dynamics
% 
% M.J. Mollema (adapted from C.C. de Visser, Delft) - 04.09.2018

    n = size(x, 1);
    xdot = zeros(n, 1);
    
    % system dynamics go here!
    xdot(1) = u(1);
    xdot(2) = u(2);
    xdot(3) = u(3);
    xdot(4) = 0;

end