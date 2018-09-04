function rankObs = kf_calcNonlinObsRank(f, h, X, X0)
% KF_CALCNONLINOBSRANK Calculates the state observability rank of the
% nonlinear given system.
%
% Inputs:
%  - f:ymbolic nonlinear state transition model 
%  - h: symbolic nonlinear observation model
%  - X: symbolic state vector
%  - X0: numeric initial state vector
%
% Output:
%  - r: rank of the observability matrix
%
% M.J. Mollema (adapted from C.C. de visser, Delft) - 04.09.2018

nstates = length(X);
nobs    = length(h);
%F = eye(size(Fx));
%Rank = [];

Hx = simplify(jacobian(h, X));
ObsMat = zeros(nobs*nstates, nstates);
Obs = sym(ObsMat, 'r');
Obs(1:nobs, :) = Hx;
Obsnum = subs(Obs, X, X0);
rankObs = double(rank(Obsnum));
fprintf('\nRank of Initial Observability matrix is %d\n', rankObs);
if (rankObs >= nstates)
    fprintf('Observability matrix is of Full Rank: the state is Observable!\n');
    return;
end
%Obs(1,:) = Hx';
LfHx = simplify(Hx * f);
for i = 2:nstates
    tic;
    LfHx = jacobian(LfHx, X);
    Obs((i-1)*nobs+1:i*nobs,:) = LfHx;
    Obsnum = subs(Obs, X, X0);
    rankObs = double(rank(Obsnum));
    fprintf('\t-> Rank of Observability matrix is %d\n', rankObs);
    if (rankObs >= nstates)
        fprintf('Observability matrix is of Full Rank: the state is Observable!\n');
        return;
    end
    LfHx = (LfHx * f);
    time = toc;
    fprintf('Loop %d took %2.2f seconds to complete\n', i, toc);
end
fprintf('WARNING: Rank of Observability matrix is %d: the state is NOT OBSERVABLE!\n', rankObs);

% fname = sprintf('NonLinObsMat_n%d', nstates);
% save(fname, 'Obs');

%Rank = [ Rank; H*F ];
%r    = rank(Rank);

end

