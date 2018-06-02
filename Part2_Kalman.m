%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2 of the AE4320:System Identification assignment
% A Kalman filter is implemented to estimate a bias in AoA
% measurements in the given dataset
%
%   Author: M.J. Mollema (adapted from original by C.C. de Visser, Delft
%   University of Technology, 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

%% Load data
filename = 'data/F16traindata_CMabV_2018';
load(filename, 'Cm', 'Z_k', 'U_k');

% transpose
Cm = Cm'; Z_k = Z_k'; U_k = U_k';

%% Set simulation parameters
dt              = 0.01;
N               = size(U_k, 2);
epsilon         = 1e-10;
doIEKF          = 1;
maxIterations   = 100;

%% Set initial values for states and statistics
Ex_0    = [Z_k(3, 1); 0; 0; 1]; % Set initial u to Vtot
x_0     = [90; -1; 1; 2];
m       = length(Ex_0);

% Initial estimate for covariance matrix
stdx_0  = [10, 10, 10, 10];
P_0     = diag(stdx_0.^2);

% System noise statistics
Ew      = [0, 0, 0, 0];
stdw    = [1e-3, 1e-3, 1e-3, 0];
Q       = diag(stdw.^2);
n       = length(stdw);
w_k     = diag(stdw) * randn(n, N) + diag(Ew) * ones(n, N);

% Measurement noise statistics
Ev      = [0, 0, 0];
stdv    = [0.01, 0.0058, 0.112];
R       = diag(stdv.^2);
nm      = length(stdv);
v_k     = diag(stdv) * randn(nm, N) + diag(Ev) * ones(nm, N);

G = eye(n);
B = eye(m);

%% Calculate batch with measurement data
% Real simulated state-variable and measurements data
% x   = x_0;
% X_k = zeros(n, N);
% Z_k = zeros(nm, N);
% U_k = zeros(m, N);

% for i = 1:N
%     dx  = kf_calc_f(0, x, U_k(:,i));
%     x   = x + (dx + w_k(:,i)) * dt;
%     
%     X_k(:,i) = x;
%     Z_k(:,i) = kf_calc_h(0, x, U_k(:,i)) + v_k(:,i);
% end

XX_k1k1     = zeros(n, N);
PP_k1k1     = zeros(n, N);
STDx_cor    = zeros(n, N);
z_pred      = zeros(nm, N);
IEKFitcount = zeros(N, 1);

x_k_1k_1 = Ex_0;
P_k_1k_1 = P_0;

%% Run EKF (with option for IEKF)
% Extended Kalman Filter (EKF)
ti = 0;
tf = dt;

printfigs = 0;
% Run the filter for all N samples
for k = 1:N
    % Prediction x(k+1|k) (integrate state eq. with one-step earlier best
    % estimate)
    [t, x_kk_1] = rk4(@kf_calc_f, x_k_1k_1, U_k(:,k), [ti, tf]);
    
    % Predicted output z(k+1|k)
    z_kk_1      = kf_calc_h(0, x_kk_1, U_k(:,k));
    z_pred(:,k) = z_kk_1;
    
    % Calculate Phi(k+1|k) and Gamma(k+1|k) (discretization)
    Fx          = kf_calc_Fx(0, x_kk_1, U_k(:,k)); % Jacobian of f(x,u)
%     [~, Psi]    = c2d(Fx, B, dt);
    [Phi, Gamma]= c2d(Fx, G, dt);
    
    % Covariance matrix of prediction P(k+1|k)
    P_kk_1      = Phi * P_k_1k_1 * Phi' + Gamma * Q * Gamma';
    P_pred      = diag(P_kk_1);
    stdx_pred   = sqrt(diag(P_kk_1));
    
    % Run IEKF if desired
    if (doIEKF)
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
            Hx  = kf_calc_Hx(0, eta1, U_k(:,k));
            
            % Check observability of the states
            if (k == 1 && itts == 1)
                check_observability
            end
            
            % The innovation matrix
            Ve  = (Hx * P_kk_1 * Hx' + R);

            % calculate the Kalman gain matrix
            K       = P_kk_1 * Hx' / Ve;
            % new observation state
            z_p     = kf_calc_h(0, eta1, U_k(:,k)) ;

            eta2    = x_kk_1 + K * (Z_k(:,k) - z_p - Hx*(x_kk_1 - eta1));
            err     = norm((eta2 - eta1), inf) / norm(eta1, inf);
        end

        IEKFitcount(k)    = itts;
        x_k_1k_1          = eta2;

    else
        % Correction
        Hx = kf_calc_Hx(0, x_kk_1, U_k(:,k)); % perturbation of h(x,u,t)
        % Pz(k+1|k) (covariance matrix of innovation)
        Ve = (Hx*P_kk_1 * Hx' + R); 

        % K(k+1) (gain)
        K = P_kk_1 * Hx' / Ve;
        % Calculate optimal state x(k+1|k+1) 
        x_k_1k_1 = x_kk_1 + K * (Z_k(:,k) - z_kk_1); 

    end    
    
    P_k_1k_1 = (eye(n) - K*Hx) * P_kk_1 * (eye(n) - K*Hx)' + K*R*K';  
    P_cor = diag(P_k_1k_1);
    stdx_cor = sqrt(diag(P_k_1k_1));

    % Next step
    ti = tf; 
    tf = tf + dt;
    
    % store results
    XX_k1k1(:,k) = x_k_1k_1;
    PP_k1k1(:,k) = diag(P_k_1k_1);
    STDx_cor(:,k) = stdx_cor;
end

% Correct alpha for bias using estimate bias
alpha = z_pred(1, :);
alpha_corr = alpha ./ (1 + XX_k1k1(4, :));

%% Plotting
F16_PlotData
