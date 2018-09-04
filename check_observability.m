function check_observabilty()
% CHECK_OBSERVABILITY Script to check if the non-linear system is
% observable using kf_calcNonlinObsRank(f, h, X, X0), which uses symbolic
% vectors as input
% 
% M.J. Mollema - 03.09.2018

%% Define symbols
syms('u', 'v', 'w','C_alpha_up', 'u_dot', 'v_dot', 'w_dot')

% State vector: x(t) = [u v w C_alpha_up]'
x  = [u; v; w; C_alpha_up];
x_0 = [100; 10; 10; 1];  % some random values

N_states = length(x);

% State transition equation: x_dot(t) = f(x(t), u(t), t) = [u(t) 0]'
f = [u_dot; v_dot; w_dot; 0];

% Measurement equation: z_n(t) = h(x(t), u(t), t) derived from definitions
% of measured/true angles and velocities
h = [atan(w/u) * (1 + C_alpha_up);
     atan(v/sqrt(u^2 + w^2));
     sqrt(u^2 + v^2 + w^2)];

% Compute rank of the observability matrix
rank = kf_calcNonlinObsRank(f, h, x, x_0);

end
