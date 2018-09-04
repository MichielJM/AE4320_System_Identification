function [XX_k1k1, z_pred, IEKFitcount] = IEKF(U_k, Z_k, dt, stdv, stdw)
% IEKF Function for running Iterated Extended Kalman Filter.
% 
% Inputs:
%  - U_k: input vector
%  - Z_k: measurement vector
%  - dt: timestep
%  - stdv: sensor noise statistics
%  - stdw: process noise statistics
% 
% Output:
%  - XX_k1k1: one-step-ahead optimal state estimation
% 
% M.J. Mollema (adapted from C.C. de Visser, Delft) - 04.09.2018

N               = size(U_k, 2);
epsilon         = 1e-10;
maxIterations   = 100;

%% Set initial values for states and statistics
Ex_0    = [Z_k(3, 1); 0.5; 0.5; 0.5]; % Set initial u to Vtot
m       = length(Ex_0);

% Initial estimate for covariance matrix
stdx_0  = [sqrt(0.1), sqrt(0.1), sqrt(0.1), sqrt(0.1)];
P_0     = diag(stdx_0.^2);

% System noise statistics
Ew      = [0, 0, 0, 0];
Q       = diag(stdw.^2);
n       = length(stdw);
w_k     = diag(stdw) * randn(n, N) + diag(Ew) * ones(n, N);

% Measurement noise statistics
Ev      = [0, 0, 0];
R       = diag(stdv.^2);
nm      = length(stdv);
v_k     = diag(stdv) * randn(nm, N) + diag(Ev) * ones(nm, N);

G = eye(n);
B = eye(m);

% Initialise array sizes
XX_k1k1     = zeros(n, N);
PP_k1k1     = zeros(n, N);
STDx_cor    = zeros(n, N);
z_pred      = zeros(nm, N);
IEKFitcount = zeros(N, 1);

% Set first estimates
x_k_1k_1 = Ex_0;
P_k_1k_1 = P_0;

ti = 0;
tf = dt;
%% Run the filter for all N samples
for k = 1:N
    % Prediction x(k+1|k) (integrate state eq. with one-step earlier best
    % estimate)
    [~, x_kk_1] = rk4(@kf_calc_f, x_k_1k_1, U_k(:,k), [ti, tf]);
    
    % Predicted output z(k+1|k)
    z_kk_1      = kf_calc_h(0, x_k_1k_1, U_k(:,k));
    z_pred(:,k) = z_kk_1;
    
    % Calculate Phi(k+1|k) and Gamma(k+1|k) (discretization)
    Fx          = kf_calc_Fx(0, x_kk_1, U_k(:,k)); % Jacobian of f(x,u)
    [Phi, Gamma]= c2d(Fx, G, dt);
    
    % Covariance matrix of prediction P(k+1|k)
    P_kk_1      = Phi * P_k_1k_1 * Phi' + Gamma * Q * Gamma';
    
    % Iterative part
    eta2    = x_kk_1;
    err     = 2 * epsilon;

    itts    = 0;
    while (err > epsilon)
        if (itts >= maxIterations)
            fprintf('Terminating IEKF: exceeded max iterations (%d)\n',...
            maxIterations);
            break
        end
        itts = itts + 1;
        eta1 = eta2;

        % Jacobian of h(x,u)
        Hx  = kf_calc_Hx(eta1, U_k(:,k));

        % The innovation matrix, P_Z(k+1|k)
        Ve  = (Hx * P_kk_1 * Hx' + R);

        % calculate the Kalman gain matrix
        K       = P_kk_1 * Hx' / Ve;
        % new observation state
        z_p     = kf_calc_h(0, eta1, U_k(:,k)) ;

        eta2    = x_kk_1 + K * (Z_k(:,k) - z_p - Hx*(x_kk_1 - eta1));
        err     = norm((eta2 - eta1), inf) / norm(eta1, inf);

        IEKFitcount(k)    = itts;
        x_k_1k_1          = eta2;
    end    
    
    P_k_1k_1 = (eye(n) - K*Hx) * P_kk_1 * (eye(n) - K*Hx)' + K*R*K'; 
    stdx_cor = sqrt(diag(P_k_1k_1));

    % Next step
    ti = tf; 
    tf = tf + dt;
    
    % store results
    XX_k1k1(:,k) = x_k_1k_1;
    PP_k1k1(:,k) = diag(P_k_1k_1);
    STDx_cor(:,k) = stdx_cor;
end

end