function test_pattern_analysis()
%TEST_PATTERN_ANALYSIS Deterministic checks for real symmetric matrices.

    % Symmetric block-tridiagonal matrix with 4x4 blocks.
    T = [1 0 0 2; 0 1 -2 0; 0 2 3 0; -2 0 0 3];
    D = diag([4 4 5 5]);
    A = sparse([D T zeros(4); T.' D T; zeros(4) T.' D]);

    bw = analyze_bandwidth(A, [0 1 2 3 4 7 8]);
    assert(bw.halfBandwidth == 7);
    assert(bw.candidateTable.OutsideNonzeros(end) == 0);
    [B,E] = split_near_band(A, 4);
    assert(isequal(A, B+E));
    assert(all(abs(nonzeros(E)) > 0));

    r4 = analyze_blocks(A, 4);
    assert(r4.isBlockTridiagonal);
    assert(r4.numNonzeroBlocks == 7);        % full 3x3 block grid
    assert(r4.numAnalyzedNonzeroBlocks == 5); % 3 diagonal + 2 above
    assert(r4.numAnalyzedBlocks == 6);
    assert(r4.blockHalfBandwidth == 1);
    assert(r4.numBlocksInRepeatedBlocks == 5); % 3 D and 2 T
    assert(r4.numTransposePairs == 3); % 3 distinct symmetric D blocks
    assert(~isfield(r4, 'numMirrorTransposePairs'));
    assert(all(arrayfun(@(g) all(g.positions(:,1) <= g.positions(:,2)), ...
        r4.equalGroups)));

    [comp, runs] = compare_block_sizes(A, [2 4 6]);
    assert(height(comp) == 3 && numel(runs) == 3);
    assert(comp.BlockTridiagonal(2));
    assert(comp.UpperNonzeroBlocks(2) == 5);
    assert(~ismember('NonMirrorTransposePairs', comp.Properties.VariableNames));

    % A nontrivial transpose relationship exists entirely above diagonal:
    % K_12 = C and K_34 = C.'; K_21 and K_43 are never analyzed.
    C = [1 2; 3 4];
    Q = sparse(8,8);
    Q(1:2,3:4) = C; Q(3:4,1:2) = C.';
    Q(5:6,7:8) = C.'; Q(7:8,5:6) = C;
    rq = analyze_blocks(Q, 2);
    assert(rq.numNonzeroBlocks == 4);
    assert(rq.numAnalyzedNonzeroBlocks == 2);
    assert(rq.numTransposePairs == 1);
    assert(numel(rq.transposeRelationships) == 1);
    pair = rq.transposeRelationships(1);
    pos = [pair.first;pair.second];
    assert(isequal(sortrows(pos), [1 2; 3 4]));

    % Partial last blocks: dimensions are part of signatures.
    partial = analyze_blocks(speye(5), 2);
    assert(partial.numNonzeroBlocks == 3);
    assert(partial.numDistinctBlocks == 2); % two 2x2 I and one 1x1

    % All-zero symmetric matrix should not cause indexing errors.
    Z = sparse(5,5);
    rz = analyze_blocks(Z, 2);
    bz = analyze_bandwidth(Z);
    assert(rz.numNonzeroBlocks == 0 && rz.isBlockTridiagonal);
    assert(rz.numAnalyzedNonzeroBlocks == 0);
    assert(rz.numTransposePairs == 0);
    assert(bz.halfBandwidth == 0 && bz.numNonzeros == 0);

    % Reject matrices that are not real symmetric.
    rejected = false;
    try
        analyze_blocks(sparse([1 1; 0 1]), 1);
    catch ME
        rejected = strcmp(ME.identifier, 'analyze_blocks:NotRealSymmetric');
    end
    assert(rejected);

    fprintf('All symmetric-matrix pattern analysis tests passed.\n');
end
