function result = analyze_blocks(A, b)
%ANALYZE_BLOCKS Statistics for contiguous b-by-b blocks of a sparse matrix.
%   R = ANALYZE_BLOCKS(A, b) reports block occupancy, block half-bandwidth,
%   repeated exact blocks, repeated support patterns, and transposed blocks.
%   A must be REAL SYMMETRIC. Repetition/transpose comparisons use only
%   the upper block triangle, including diagonal blocks (p <= q).
%   Mirrored pairs (p,q)/(q,p) are NEVER compared or reported.
%
%   The last block row/column may be smaller than b. Two blocks are counted
%   equal only when their actual dimensions agree. Equality is EXACT (not
%   tolerant): useful for first-pass structural discovery on supplied data.
%   All nonzero blocks, not all possible blocks, are processed. No full(A).
%
%   A group contains the positions [blockRow, blockColumn] of its members.
%   transposeRelationships describes equal-transpose groups and counts pairs
%   of DISTINCT upper-triangle blocks. Full-matrix occupancy counts are
%   retained independently for bandwidth and block-tridiagonal testing.

    if ~ismatrix(A) || ~isnumeric(A)
        error('analyze_blocks:InvalidMatrix', 'A must be a numeric matrix.');
    end
    validateattributes(b, {'numeric'}, ...
        {'scalar','real','finite','positive','integer'}, mfilename, 'b');
    b = double(b);
    [m, n] = size(A);
    if m ~= n || ~isreal(A) || ~issymmetric(A)
        error('analyze_blocks:NotRealSymmetric', ...
            'This version requires a real symmetric square matrix A.');
    end
    nBR = ceil(m / b);
    nBC = ceil(n / b);

    [i, j, v] = find(A);
    br = floor((double(i) - 1) / b) + 1;
    bc = floor((double(j) - 1) / b) + 1;
    localR = mod(double(i) - 1, b) + 1;
    localC = mod(double(j) - 1, b) + 1;

    blockNNZ = sparse(br, bc, ones(numel(i),1), nBR, nBC);
    blockOccupancy = spones(blockNNZ);
    [allOccupiedR, allOccupiedC] = find(blockOccupancy);
    numberOfOccupiedAll = numel(allOccupiedR);

    % Only upper-triangle blocks participate in pattern comparisons.
    [occupiedR, occupiedC] = find(triu(blockOccupancy));
    numberOfOccupied = numel(occupiedR);
    keep = br <= bc;
    br = br(keep);
    bc = bc(keep);
    localR = localR(keep);
    localC = localC(keep);
    v = v(keep);

    % Group upper-triangle nonzeros by block, in column-major order.
    if isempty(v)
        members = cell(0,1);
    else
        linearBlockId = sub2ind([nBR, nBC], br, bc);
        [~, ~, blockId] = unique(linearBlockId);
        members = accumarray(blockId, (1:numel(v)).', ...
            [numberOfOccupied, 1], @(x){x}, {[]});
    end

    % Keys are exact and collision-free for double-valued data. Values use
    % IEEE-754 bit representations (real and imaginary components).
    supportMap = containers.Map('KeyType','char','ValueType','any');
    exactMap = containers.Map('KeyType','char','ValueType','any');
    records = repmat(struct('row',0,'col',0,'height',0,'width',0, ...
        'nonzeros',0,'supportKey','','exactKey','','transposeKey',''), ...
        numberOfOccupied, 1);

    for k = 1:numberOfOccupied
        idx = members{k};
        p = occupiedR(k);
        q = occupiedC(k);
        nr = min(b, m - (p-1)*b);
        nc = min(b, n - (q-1)*b);
        locations = [localR(idx), localC(idx)];
        [locations, order] = sortrows(locations, [1,2]);
        values = v(idx(order));

        supportKey = makeSupportKey(nr, nc, locations);
        exactKey = makeValueKey(nr, nc, locations, values);
        transposeKey = makeValueKey(nc, nr, locations(:,[2,1]), values);

        records(k).row = p;
        records(k).col = q;
        records(k).height = nr;
        records(k).width = nc;
        records(k).nonzeros = numel(idx);
        records(k).supportKey = supportKey;
        records(k).exactKey = exactKey;
        records(k).transposeKey = transposeKey;
        supportMap = appendIndex(supportMap, supportKey, k);
        exactMap = appendIndex(exactMap, exactKey, k);
    end

    [supportGroups, repeatedSupportCount, supportMembers] = ...
        groupedRepetitions(supportMap, records);
    [equalGroups, repeatedEqualCount, equalMembers] = ...
        groupedRepetitions(exactMap, records);

    % Count transpose pairs at the signature-group level rather than
    % materializing every pair (which may grow quadratically).
    relationshipTemplate = struct('first',zeros(0,2), ...
        'second',zeros(0,2),'numberOfPairs',0);
    relationships = repmat(relationshipTemplate, 0, 1);
    visited = containers.Map('KeyType','char','ValueType','logical');
    exactKeys = keys(exactMap);
    numPairs = 0;
    for t = 1:numel(exactKeys)
        key = exactKeys{t};
        if isKey(visited, key)
            continue;
        end
        visited(key) = true;
        group1 = exactMap(key);
        other = records(group1(1)).transposeKey;
        if ~isKey(exactMap, other)
            continue;
        end
        group2 = exactMap(other);
        visited(other) = true;
        if strcmp(key, other)
            pairCount = numel(group1) * (numel(group1) - 1) / 2;
        else
            pairCount = numel(group1) * numel(group2);
        end
        if pairCount == 0
            continue;
        end
        numPairs = numPairs + pairCount;
        r = relationshipTemplate;
        r.first = [[records(group1).row].', [records(group1).col].'];
        r.second = [[records(group2).row].', [records(group2).col].'];
        r.numberOfPairs = pairCount;
        relationships(end+1,1) = r; %#ok<AGROW>
    end

    result = struct();
    result.blockSize = b;
    result.matrixSize = [m n];
    result.blockGridSize = [nBR nBC];
    result.numBlocks = nBR * nBC;
    result.numNonzeroBlocks = numberOfOccupiedAll; % whole matrix
    result.numAnalyzedBlocks = nBR * (nBR + 1) / 2;
    result.numAnalyzedNonzeroBlocks = numberOfOccupied; % upper triangle
    result.blockOccupancy = blockOccupancy;
    result.blockNonzeroCounts = blockNNZ;
    result.blockHalfBandwidth = 0;
    if numberOfOccupiedAll > 0
        result.blockHalfBandwidth = max(abs(allOccupiedR - allOccupiedC));
    end
    result.isBlockTridiagonal = result.blockHalfBandwidth <= 1;
    result.numDistinctSupports = supportMap.Count;
    result.numRepeatedSupportGroups = repeatedSupportCount;
    result.numBlocksInRepeatedSupports = supportMembers;
    result.numDistinctBlocks = exactMap.Count;
    result.numRepeatedBlockGroups = repeatedEqualCount;
    result.numBlocksInRepeatedBlocks = equalMembers;
    result.supportGroups = supportGroups;
    result.equalGroups = equalGroups;
    result.transposeRelationships = relationships;
    result.numTransposePairs = numPairs;
end

function key = makeSupportKey(nr, nc, locations)
    locations = sortrows(locations, [1,2]);
    key = sprintf('%dx%d|', nr, nc);
    if ~isempty(locations)
        key = [key, sprintf('%d,%d;', locations.')]; %#ok<AGROW>
    end
end

function key = makeValueKey(nr, nc, locations, values)
    [locations, order] = sortrows(locations, [1,2]);
    values = double(values(order));
    key = makeSupportKey(nr, nc, locations);
    realHex = reshape(num2hex(real(values)).', 1, []);
    imagHex = reshape(num2hex(imag(values)).', 1, []);
    key = [key 'V|' realHex '|' imagHex];
end

function dict = appendIndex(dict, key, index)
    if isKey(dict,key)
        dict(key) = [dict(key), index];
    else
        dict(key) = index;
    end
end

function [groups, nRepeatGroups, membersOfRepeatGroups] = groupedRepetitions(dict, records)
    groupTemplate = struct('count',0,'positions',zeros(0,2));
    groups = repmat(groupTemplate, 0, 1);
    nRepeatGroups = 0;
    membersOfRepeatGroups = 0;
    allKeys = keys(dict);
    for k = 1:numel(allKeys)
        members = dict(allKeys{k});
        if numel(members) < 2
            continue;
        end
        group = groupTemplate;
        group.count = numel(members);
        group.positions = [[records(members).row].', ...
                           [records(members).col].'];
        groups(end+1,1) = group; %#ok<AGROW>
        nRepeatGroups = nRepeatGroups + 1;
        membersOfRepeatGroups = membersOfRepeatGroups + numel(members);
    end
    if ~isempty(groups)
        [~, order] = sort([groups.count], 'descend');
        groups = groups(order);
    end
end
