function write_block_report(report, outputPath, maxGroups)
%WRITE_BLOCK_REPORT Human-readable summary of repeated block relationships.
    if nargin < 3, maxGroups = 15; end
    fid = fopen(outputPath, 'w');
    if fid < 0
        error('write_block_report:CannotOpen', 'Cannot open %s', outputPath);
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, 'BLOCK SIZE: %d\n', report.blockSize);
    fprintf(fid, 'Block grid: %d x %d\n', report.blockGridSize);
    fprintf(fid, 'Nonzero blocks (whole matrix): %d / %d\n', ...
        report.numNonzeroBlocks, report.numBlocks);
    fprintf(fid, 'Analyzed blocks (upper triangle): %d / %d\n', ...
        report.numAnalyzedNonzeroBlocks, report.numAnalyzedBlocks);
    fprintf(fid, 'Block half-bandwidth: %d\n', report.blockHalfBandwidth);
    fprintf(fid, 'Block tridiagonal: %d\n', report.isBlockTridiagonal);
    fprintf(fid, 'Distinct supports: %d\n', report.numDistinctSupports);
    fprintf(fid, 'Distinct exact blocks: %d\n', report.numDistinctBlocks);
    fprintf(fid, 'Transpose pairs (distinct upper-triangle blocks only): %d\n', ...
        report.numTransposePairs);

    fprintf(fid, '\nREPEATED SUPPORT PATTERNS\n');
    printGroups(fid, report.supportGroups, maxGroups);
    fprintf(fid, '\nEXACTLY EQUAL BLOCKS\n');
    printGroups(fid, report.equalGroups, maxGroups);
    fprintf(fid, '\nTRANSPOSE RELATIONSHIPS (distinct upper-triangle blocks only)\n');
    pairs = report.transposeRelationships;
    if isempty(pairs)
        fprintf(fid, '(none)\n');
    end
    [~, order] = sort([pairs.numberOfPairs], 'descend');
    pairs = pairs(order);
    for k = 1:min(maxGroups, numel(pairs))
        fprintf(fid, '%d distinct block pairs: %s  <->  %s\n', ...
            pairs(k).numberOfPairs, ...
            positionsText(pairs(k).first, 8), ...
            positionsText(pairs(k).second, 8));
    end
    if numel(pairs) > maxGroups
        fprintf(fid, '... %d additional relationships omitted\n', ...
            numel(pairs)-maxGroups);
    end
end

function printGroups(fid, groups, maxGroups)
    if isempty(groups)
        fprintf(fid, '(none)\n');
        return;
    end
    for k = 1:min(maxGroups, numel(groups))
        fprintf(fid, '%d blocks: %s\n', groups(k).count, ...
            positionsText(groups(k).positions, 15));
    end
    if numel(groups) > maxGroups
        fprintf(fid, '... %d additional groups omitted\n', ...
            numel(groups)-maxGroups);
    end
end

function txt = positionsText(pos, maxCount)
    n = min(size(pos,1), maxCount);
    pieces = cell(n,1);
    for k = 1:n
        pieces{k} = sprintf('(%d,%d)', pos(k,1), pos(k,2));
    end
    txt = strjoin(pieces, ' ');
    if size(pos,1) > maxCount
        txt = [txt sprintf(' ... (+%d)',size(pos,1)-maxCount)];
    end
end
