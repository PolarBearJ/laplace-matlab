% July 20th, 2022 Joshua Stevens 20907520
% laplace2d
% The approximation of the solution to Laplace’s equation will give the
% steady state of the function u, that is, the center of the boundary will
% be equal to the average of the boundaries.
%
% Parameters
% ==========
%    U - The boundary matrix in two or three dimensions.
%
% Return Values
% =============
%    U_out - The steady state matrix of the boundary, and solution the Laplace equation.

function [U_out] = laplace2d( U )
% Error and Warning Checking
% ==========================
%   
%   Check if U is a matrix, then
%   ensure the dimensions of the boundary matrix. If we are solving in 2
%   or 3 dimensions, ensure the dimensions of the boundary matrix match.
    if ~ismatrix( U )
        throw( MException( 'MATLAB:invalid_argument', ...
            'the argument U is not a matrix' ) );
    end
    if ~ndims( U ) == 2
            throw( MException( 'MATLAB:invalid_argument', ...
            'the argument U is not a matrix of 2 demensions,  U is of %d dimensions', ndims(U) ) );
    end
% Initialization
% ==============
%
%   Initialize number of x and y coordinates using the size of the matrix
%   given in each dimension.
    [n_x, n_y] = size( U );
    U_out = U;
% Mapping the unknown points to a unique number from 1 to m
% =========================================================
%
%   First create a counter variable m set to zero, then create two matrices,
%   one of which (u to w) will take a x and y coordinate (and z if in 3D)
%   and map it to its corresponding w value, the matrix should be of size
%   x coordinates by y coordinates (again, by z coordinates as well if in 3D).
%   The other (w to u) will take a w value and map it to its corresponding
%   coordinates, this matrix should be of size 2 by number of coordinates
%   in each dimension multiplied by each other. This is done by looping
%   through each coordinate value of the boundary matrix, and at each value
%   check if this value is equal to negative infinity. If this condition is
%   true, we increment the m value, set the value of u to w at this point
%   to m and set the value of w to u at m to the coordinates.
    u_to_w = zeros( n_x, n_y );
    w_to_u = zeros( 2, n_x * n_y );
    m = 0;
    for ix = 1:n_x
        for iy = 1:n_y
            if U(ix, iy) == -Inf
                m = m + 1;
                u_to_w(ix, iy) = m;
                w_to_u(:, m) = [ix, iy]';
            end
        end
    end
% Creating and solving a system of linear equations
% =================================================
%
%   Create a zero matrix M m by m in size and a zero
%   vector b of length m. Create a loop that iterates through each w value,
%   for each check if the following conditions in each 4 (or 6 if 3D)
%   adjacent are met.
    M = zeros( m, m );
    b = zeros( m, 1 );

    for i = 1:m
        c = w_to_u(:, i);
        c_adjacent = [ c(1) - 1, c(2) + 0
                       c(1) + 1, c(2) + 0
                       c(1) + 0, c(2) - 1
                       c(1) + 0, c(2) + 1];
        for j = 1:4
            % Find the value of the adjacent
            adjacent_value = U(c_adjacent(j,1), c_adjacent(j,2));
            % If Dirichlet boundary condition
            % Subtract 1 from the ith diagonal entry of M (current diagonal position)
            % Subtract the value from the ith entry of the vector b 
            % (subtract the value given by the boundary at the position of the vector b)
            if ~isnan(adjacent_value) && adjacent_value ~= -Inf
                M(i, i) = M(i, i) - 1;
                b(i) = b(i) - adjacent_value;
            % A different unknown (-Inf) point in this case, the wth unknown
            % Subtract 1 from the ith diagonal entry of M
            % Add 1 to the (i, w)th entry of M
            elseif ~isnan(adjacent_value)
                M(i, i) = M(i, i) - 1;
                M(i, u_to_w(c_adjacent(j, 1), c_adjacent(j,2))) = M(i, u_to_w(c_adjacent(j, 1), c_adjacent(j,2))) + 1;
            end
        end
    end
% Substituting the values back into the matrix U_out
% ===================================================
%
%   Cycle through each solution value, find the coordinate of each and
%   insert the approximated value in the coordinate position. Coordinates
%   are found using the w to u function.

    % Once matrix M and vector b are filled, solve for u
    u = M \ b;
    % Cycle through each solution value, find the coordinate of each
    % and insert the approximated value in the correct position
    for k = 1:m
        coords = w_to_u(:, k);
        U(coords(1), coords(2)) = u(k);
    end    
    U_out = U;
end
