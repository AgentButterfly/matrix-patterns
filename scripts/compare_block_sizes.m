function [comparison, analyses] = compare_block_sizes(A, blockSizes)
%COMPARE_BLOCK_SIZES Compare structural patterns for several block sizes.
%   [TABLE, RESULTS] = COMPARE_BLOCK_SIZES(A, [2 4 6 8]);
%   A must be real symmetric. Occupancy metrics use the full block grid;
%   repetitions and transpose comparisons use ONLY upper-triangle blocks.
%   The table is for manual exploration, not a circuit-cost ranking.

    validateattributes(blockSizes, {'numeric'}, ...
        {'vector','real','finite','positive','integer'}, ...
        mfilename, 'blockSizes');
    blockSizes = unique(double(blockSizes(:)), 'stable');
    if isempty(blockSizes)
        error('compare_block_sizes:NoSizes', 'Specify at least one block size.');
    end

    k = numel(blockSizes);
    analyses = cell(k, 1);
    nonzeroBlocks = zeros(k,1);
    totalBlocks = zeros(k,1);
    upperNonzeroBlocks = zeros(k,1);
    upperTotalBlocks = zeros(k,1);
    blockHalfBandwidths = zeros(k,1);
    blockTridiagonal = false(k,1);
    distinctSupports = zeros(k,1);
    blocksWithRepeatedSupports = zeros(k,1);
    distinctExactBlocks = zeros(k,1);
    blocksWithExactRepeats = zeros(k,1);
    transposePairs = zeros(k,1);

    for t = 1:k
        r = analyze_blocks(A, blockSizes(t));
        analyses{t} = r;
        nonzeroBlocks(t) = r.numNonzeroBlocks;
        totalBlocks(t) = r.numBlocks;
        upperNonzeroBlocks(t) = r.numAnalyzedNonzeroBlocks;
        upperTotalBlocks(t) = r.numAnalyzedBlocks;
        blockHalfBandwidths(t) = r.blockHalfBandwidth;
        blockTridiagonal(t) = r.isBlockTridiagonal;
        distinctSupports(t) = r.numDistinctSupports;
        blocksWithRepeatedSupports(t) = r.numBlocksInRepeatedSupports;
        distinctExactBlocks(t) = r.numDistinctBlocks;
        blocksWithExactRepeats(t) = r.numBlocksInRepeatedBlocks;
        transposePairs(t) = r.numTransposePairs;
    end

    comparison = table(blockSizes, nonzeroBlocks, totalBlocks, ...
        upperNonzeroBlocks, upperTotalBlocks, blockHalfBandwidths, ...
        blockTridiagonal, distinctSupports, blocksWithRepeatedSupports, ...
        distinctExactBlocks, blocksWithExactRepeats, transposePairs, ...
        'VariableNames', {'BlockSize','NonzeroBlocks','TotalBlocks', ...
        'UpperNonzeroBlocks','UpperTotalBlocks','BlockHalfBandwidth', ...
        'BlockTridiagonal','DistinctSupports','BlocksInRepeatedSupports', ...
        'DistinctExactBlocks','BlocksInExactRepeats','TransposePairs'});
end
