
% BSPLINEN_CART2BARY converts cartesian coordinates in n-space to barycentric
%   coordinates in n+1 space.
%
%              Copyright: C.C. de Visser, Delft University of Technology, 2007
%              email: c.c.devisser@tudelft.nl
%   
%   LAMBDA = BSPLINEN_CART2BARY(SIMPLEX, X) converts the cartesian coordinates of X
%   into barycentric coordinates LAMBDA with respect to SIMPLEX.
%
%   Input parameters are:
%       SIMPLEX is a matrix containing the vertex coordinates of a simplex in n space.
%           The rows in SIMPLEX are the vertices, the columns the coordinates.
%           Therefore SIMPLEX always has n+1 rows and n columns.
%       X is a maxtrix holding points in n-space in cartesian coordinates. 
%           X must have n columns, with n the dimension of SIMPLEX.
%
%   Output from BSPLINEN_CART2BARY is:
%
%       LAMBDA the barycentric coordinates of X with respect to SIMPLEX. When 
%           LAMBDA contains only positive values, X is within the convex
%           hull of SIMPLEX. If LAMBDA contains one or more negative
%           values, then X is outside of the convex hull of SIMPLEX.
function Lambda = bsplinen_cart2bary(simplex, X)
    
    % The reference vertex is always chosen as the first simplex vertex.
    % This can be done because barycentric coordinates are not dependent on
    % the reference point.
    v0      = simplex(1, :);
    vcount2 = size(simplex, 1) - 1;
    Xcount  = size(X, 1);
    
    Lambda = zeros(Xcount, vcount2 + 1);
    %vcount2 = length(simplex(:, 1)) - 1;
    
    % assemble matrix A
    A = zeros(vcount2, vcount2);
    count = 1;
    for i = 2:length(simplex(:, 1))
        A(:, count) = (simplex(i, :) - v0)';
        count = count + 1;
    end
    
    
    for i = 1:Xcount
        % relative coordinates of x
        p = (X(i, :) - v0)';

        % the last (n) barycentric coordinates. 
        lambda1 = A \ p;

        % the first barycentric coordinate; lambda0
        lambda0 = 1 - sum(lambda1);

        % insert lambda0 into the Lambda vector
        Lambda(i, 1) = lambda0;
        Lambda(i, 2:end) = lambda1;
    end




