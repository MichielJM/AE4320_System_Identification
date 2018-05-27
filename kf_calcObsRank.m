%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r = fpr_calcObsRank(H, Fx) Calculates the state observability rank
%       of the discretized linear system given as:
%           xdot = Fx*x + G*u
%           zm   = H*x + D*u
%
%   Input arguments are:
%         Fx - state transition matrix 
%         H  - observation matrix
%
%   Output argument is:
%         r - the rank of the observability matrix
%   
%   Authors: C.C. de Visser, Delft University of Technology, 2011
%            J. Oliveira, Delft University of Technology, 2004
%   email: c.c.devisser@tudelft.nl
%   Version: 2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = kf_calcObsRank(H, Fx)

    nstates = size(Fx,1);

    F = eye(size(Fx));
    Rank = [];
    for i = 1:nstates-1,
       Rank = [ Rank; H*F ];
       F = F*Fx;
    end
    Rank = [ Rank; H*F ];
    r    = rank(Rank);

