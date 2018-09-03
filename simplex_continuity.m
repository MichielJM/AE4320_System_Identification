function H = simplex_continuity(tri, spline_poly_order, spline_cont_order,...
    c_OLS_coeff)
% SIMPLEX_CONTINUITY Creates the continuity matrix H
% 
% Inputs:
%  - tri: MATLAB triangulation object
%  - spline_poly_order: order of the polynomial spline
%  - spline_cont_order: order to which the splines have to be continuous on
%       the triangle edges
%  - c_OLS_coeff: 'symbolic' b-coefficient vector
% 
% Outputs:
%  - H: continuity matrix
% 
% M.J. Mollema - 03.09.2018

%% Loop over all orders and edges to determine continuity matrix H
% Find all edges
triangles   = sort(tri.ConnectivityList, 2);
edgess      = edges(tri);
boundary    = freeBoundary(tri);
int_edges   = setdiff(edgess, [boundary; fliplr(boundary)], 'rows'); % Find only interior edges

% Generate exponentials
exponentials    = gen_exp(3, spline_poly_order);

% Empty array to store data
H = [];

for i = 1:size(int_edges, 1)
    
    for order = 0:spline_cont_order
        
        % Determine adjacent triangles of current edge and Out-Of-Edge
        % vertices.
        edge = int_edges(i, :);
        current_triangles = find(sum((triangles == edge(1)) + (triangles == edge(2)), 2) == 2);
        OOE_vertex_1 = triangles(current_triangles(1), :);
        OOE_vertex_1((OOE_vertex_1 == edge(1)) | (OOE_vertex_1 == edge(2))) = [];
        OOE_vertex_2 = triangles(current_triangles(2), :);
        OOE_vertex_2((OOE_vertex_2 == edge(1)) | (OOE_vertex_2 == edge(2))) = [];
        
        % Find barycentric coords of OOE_1 wrt triangle 2
        simplex_2 = tri.Points(triangles(current_triangles(2), :), :);
        OOE_bary_12 = bsplinen_cart2bary(simplex_2, tri.Points(OOE_vertex_1, :));
        
        % Find location of single non-zero value in vertex multi-index
        single_nonzero_loc_1 = find(triangles(current_triangles(1), :) == OOE_vertex_1);
        single_nonzero_loc_2 = find(triangles(current_triangles(2), :) == OOE_vertex_2);

        % Get all permutations of exponentials which have value of 'order'
        % in their single nonzero location
        indices     = exponentials(:, single_nonzero_loc_1) == order;
        LH_part     = exponentials(indices, :);
        
        % Gamma permutations
        gamma = gen_exp(3, order);
        
        % Right hand part
        RH_part = [];
        % Loop over all LH_parts and per part add all gamma permutations
        % and create continuity matrix H
        for j = 1:size(LH_part, 1)
            
            % Get k0 and k1 from left hand part
            ks = LH_part(j, :);
            ks(:, single_nonzero_loc_1) = [];
            RH = zeros(1, 3);
            if single_nonzero_loc_2 == 1
                RH(2:3) = ks;
            elseif single_nonzero_loc_2 == 2
                RH(1) = ks(1);
                RH(3) = ks(2);
            else
                RH(1:2) = ks;
            end
            
            % All permutations of RH + gamma
            RH = RH + gamma;
            RH_part = vertcat(RH_part, RH);
            
            % Smoothness conditions (find index of corresponding entry in
            % c_OLS and set that to -1 or b(v))
            h_vector = zeros(1, size(c_OLS_coeff, 1));
            % Left side
            LH = LH_part(j, :);
            LH(single_nonzero_loc_1) = order; % Set to order at single-nonzero-location
            LH_idx = find(ismember(c_OLS_coeff, [current_triangles(1), LH], 'rows'));
            h_vector(LH_idx) = -1;
            % Right side
            RH_idx = find(ismember(c_OLS_coeff, [current_triangles(2)*ones(size(RH, 1), 1), RH], 'rows'));
            b_v = prod(OOE_bary_12.^gamma, 2);
            h_vector(RH_idx) = b_v;
            
            % Fill H
            H = vertcat(H, h_vector);
            
        end
            
    end
    
end

end