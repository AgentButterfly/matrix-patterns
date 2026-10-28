function result = analyze_bandwidth(A, candidateBandwidths)
%ANALYZE_BANDWIDTH Exact diagonal-offset statistics and near-banded splits.
%   R = ANALYZE_BANDWIDTH(A) chooses a small set of candidate half-bandwidths.
%   R = ANALYZE_BANDWIDTH(A, H) uses the nonnegative integer widths in H.
%
%   For each h, A = B_h + E_h, with B_h retaining entries |i-j| <= h.
%   candidateTable reports nnz(E_h), coverage, and the number kept in B_h.
%   Both this function and ANALYZE_BLOCKS operate on the original indexing.
%
%   No dense copy of A is created.

    if nargin < 2
        candidateBandwidths = [];
    end
    if ~ismatrix(A) || ~isnumeric(A)
        error('analyze_bandwidth:InvalidMatrix', 'A must be a numeric matrix.');
    end

    [i, j] = find(A);
    offsets = double(j) - double(i);
    distances = abs(offsets);
    numNZ = numel(distances);

    if numNZ == 0
        halfBandwidth = 0;
        offsetTable = table(zeros(0,1), zeros(0,1), ...
            'VariableNames', {'Offset','Nonzeros'});
    else
        halfBandwidth = max(distances);
        [distinctOffsets, ~, offsetIndex] = unique(offsets);
        counts = accumarray(offsetIndex, 1);
        offsetTable = table(distinctOffsets(:), counts(:), ...
            'VariableNames', {'Offset','Nonzeros'});
    end

    if isempty(candidateBandwidths)
        candidates = [0 1 2 3 4 6 8 12 16 24 32 64 halfBandwidth];
        if numNZ > 0
            sortedDistances = sort(distances);
            ranks = max(1, ceil([0.90 0.95 0.99] * numNZ));
            candidates = [candidates sortedDistances(ranks).'];
        end
        candidates = unique(min(candidates, halfBandwidth));
    else
        validateattributes(candidateBandwidths, {'numeric'}, ...
            {'vector','real','finite','nonnegative','integer'}, ...
            mfilename, 'candidateBandwidths');
        candidates = unique(double(candidateBandwidths(:).'));
    end

    candidates = candidates(:);
    outside = zeros(numel(candidates), 1);
    for k = 1:numel(candidates)
        outside(k) = nnz(distances > candidates(k));
    end
    inside = numNZ - outside;
    if numNZ == 0
        coverage = ones(size(inside));
    else
        coverage = inside ./ numNZ;
    end

    result = struct();
    result.matrixSize = size(A);
    result.numNonzeros = numNZ;
    result.halfBandwidth = halfBandwidth;
    result.offsetTable = offsetTable;
    result.candidateTable = table(candidates, inside, outside, coverage, ...
        'VariableNames', {'HalfBandwidth','InsideNonzeros', ...
                          'OutsideNonzeros','Coverage'});
end
