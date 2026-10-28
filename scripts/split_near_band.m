function [B, E] = split_near_band(A, h)
%SPLIT_NEAR_BAND Exact sparse decomposition into banded core and correction.
%   [B, E] = SPLIT_NEAR_BAND(A, h) gives A = B + E, with B(i,j) = A(i,j)
%   for |i-j| <= h, and E containing only entries outside this band.
    validateattributes(h, {'numeric'}, ...
        {'scalar','real','finite','nonnegative','integer'}, mfilename, 'h');
    if ~ismatrix(A) || ~isnumeric(A)
        error('split_near_band:InvalidMatrix', 'A must be numeric and 2-D.');
    end
    [m,n] = size(A);
    [i,j,v] = find(A);
    inside = abs(double(i)-double(j)) <= h;
    B = sparse(i(inside), j(inside), v(inside), m, n);
    E = sparse(i(~inside), j(~inside), v(~inside), m, n);
end
