


% BSPLINEN_BARY2CART converts barycentric coordinates in n+1 spacec to
%   cartesian coordinates in n-space.
%   
%   X = BSPLINEN_BARY2CART(SIMPLEX, LAMBDA) converts the barycentric coordinates of
%   LAMBDA with respect to the vertices of SIMPLEX to (global)
%   cartesian coordinates X
%
%   Input parameters are:
%       SIMPLEX contains the vertex coordinates of a simplex in n space.
%           The rows in SIMPLEX are the vertices, the columns the coordinates.
%           Therefore SIMPLEX always has n+1 rows and n columns.
%
%       LAMBDA a matrix holding barycentric coordinates with respect to SIMPLEX. 
%           LAMBDA must have n+1 columns, with n the dimension of SIMPLEX.
%
%   Output from BSPLINEN_BARY2CART is:
%
%       X the global cartesian coordinates of the point LAMBDA. When 
%           LAMBDA contains only positive values, X is within the convex
%           hull of SIMPLEX. If LAMBDA contains one or more negative
%           values, then X is outside of the convex hull of SIMPLEX.
function X = bsplinen_bary2cart(simplex, lambda)

    Lcount = size(lambda, 1);
    X = zeros(Lcount, size(lambda, 2) - 1);
    
    v0 = simplex(1, :);
    vcount2 = size(simplex, 1);
    
    for j = 1:Lcount
        for i = 2:vcount2
            X(j, :) = X(j, :) + lambda(j, i) * (simplex(i, :) - v0);
        end
        X(j, :) = X(j, :) + v0;
    end
    
    