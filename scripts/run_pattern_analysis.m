%% Sparse-matrix pattern analysis — edit these settings
% This script is designed to live in matrix-patterns/scripts/.
group = 'HB';
name = 'bcsstk16';
blockSizes = [2 3 4 6 8 12 14 16 20 22 24 26 28 30 32 34 36 38 40 48 56 64 72 74 76 78 80 88 96 104 112 120 128 136 138 140 142 144 ];
candidateBandwidths = [];  % [] = automatic candidates; or e.g. [0 2 4 6 8 12]
bandToSplit = [];           % e.g. 4 to save A = B_4 + E_4; [] = no split
showMatrixTable = false;   % Avoid full(A) for large SuiteSparse matrices
saveBlockPlots = true;

%% Repository-relative paths
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
addpath(scriptDir);
dataDir = fullfile(repoDir, 'data', 'raw', group);
imageDir = fullfile(repoDir, 'images', group);
reportDir = fullfile(repoDir, 'results', group, name);
if ~isfolder(dataDir), mkdir(dataDir); end
if ~isfolder(imageDir), mkdir(imageDir); end
if ~isfolder(reportDir), mkdir(reportDir); end

filePath = fullfile(dataDir, [name '.mat']);
imagePath = fullfile(imageDir, [name '.png']);

%% Load the matrix (cached MAT or NIST Harwell-Boeing fallback)
% No SuiteSparse websave call: its file server is currently returning 403.
% The helper caches NIST downloads in data/raw/HB/ for reuse.
A = load_matrix_with_fallback(group, name, dataDir);
if size(A,1) ~= size(A,2) || ~isreal(A) || ~issymmetric(A)
    error('run_pattern_analysis:NotRealSymmetric', ...
        'The block detector currently requires a real symmetric matrix.');
end
fprintf('\nMatrix: %s/%s | %d x %d | nnz = %d\n', ...
    group, name, size(A,1), size(A,2), nnz(A));

%% Analyze the matrix
bandReport = analyze_bandwidth(A, candidateBandwidths);
[comparison, blockReports] = compare_block_sizes(A, blockSizes);

fprintf('\nHalf-bandwidth: %d\n', bandReport.halfBandwidth);
disp(bandReport.candidateTable);
fprintf('\nBlock-size comparison:\n');
disp(comparison);

%% Save numerical reports
save(fullfile(reportDir, 'analysis.mat'), ...
    'bandReport', 'comparison', 'blockReports', 'blockSizes', ...
    '-v7.3');
writetable(bandReport.offsetTable, ...
    fullfile(reportDir, 'diagonal_offsets.csv'));
writetable(bandReport.candidateTable, ...
    fullfile(reportDir, 'candidate_bands.csv'));
writetable(comparison, fullfile(reportDir, 'block_sizes.csv'));
for t = 1:numel(blockSizes)
    write_block_report(blockReports{t}, ...
        fullfile(reportDir, sprintf('blocks_b%d.txt', blockSizes(t))));
end

%% Optional exact near-band decomposition
if ~isempty(bandToSplit)
    [B, E] = split_near_band(A, bandToSplit);
    fprintf('Band h=%d: %d retained entries; %d exceptional entries.\n', ...
        bandToSplit, nnz(B), nnz(E));
    save(fullfile(reportDir, sprintf('split_h%d.mat',bandToSplit)), ...
        'B', 'E', '-v7.3');
    f = figure('Visible','off');
    spy(B);
    title(sprintf('%s/%s — band h=%d',group,name,bandToSplit), ...
          'Interpreter','none');
    exportgraphics(f, fullfile(reportDir,sprintf('band_h%d.png',bandToSplit)), ...
        'Resolution',200);
    close(f);
    f = figure('Visible','off');
    spy(E);
    title(sprintf('%s/%s — exceptions beyond h=%d',group,name,bandToSplit), ...
          'Interpreter','none');
    exportgraphics(f, fullfile(reportDir,sprintf('exceptions_h%d.png',bandToSplit)), ...
        'Resolution',200);
    close(f);
end

%% Save the ordinary sparsity plot
fig = figure('Visible','on');
spy(A);
title([group '/' name], 'Interpreter', 'none');
exportgraphics(fig, imagePath, 'Resolution', 300);
fprintf('Saved image: %s\n', imagePath);

%% Save block-occupancy plots for all candidate block sizes
if saveBlockPlots
    for t = 1:numel(blockSizes)
        b = blockSizes(t);
        r = blockReports{t};
        f = figure('Visible','off');
        spy(r.blockOccupancy);
        title(sprintf('%s/%s — block occupancy, b = %d', group, name, b), ...
              'Interpreter', 'none');
        xlabel('Block column'); ylabel('Block row');
        outPath = fullfile(reportDir, sprintf('blocks_b%d.png', b));
        exportgraphics(f, outPath, 'Resolution', 200);
        close(f);
    end
end

%% Optional small-matrix table (dense conversion only when reasonable)
if showMatrixTable
    if numel(A) <= 40000
        f = uifigure('Name', [group '/' name]);
        uitable(f, 'Data', full(A), 'Position', [20 20 700 450]);
        f.Position = [100 100 740 490];
    else
        warning('Matrix too large for a dense table; skipping display.');
    end
end
fprintf('Saved analysis to: %s\n', reportDir);
