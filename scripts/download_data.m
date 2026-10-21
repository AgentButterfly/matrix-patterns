
%% Choose the matrix
group = 'HB';
name  = 'bcsstk04';

%% Set paths relative to this script
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);

dataDir = fullfile(repoDir, 'data', 'raw', group);
imageDir = fullfile(repoDir, 'images', group);

if ~isfolder(dataDir)
    mkdir(dataDir);
end

if ~isfolder(imageDir)
    mkdir(imageDir);
end

filePath = fullfile(dataDir, [name '.mat']);
imagePath = fullfile(imageDir, [name '.png']);

%% Download the matrix if it is not already saved
url = sprintf('https://sparse.tamu.edu/mat/%s/%s.mat', ...
              group, name);

if ~isfile(filePath)
    fprintf('Downloading %s/%s...\n', group, name);

    options = weboptions('Timeout', 120);
    websave(filePath, url, options);
else
    fprintf('Using existing file: %s\n', filePath);
end

%% Load the matrix
loaded = load(filePath, 'Problem');
Problem = loaded.Problem;
A = Problem.A;

%% Display basic information
fprintf('Matrix: %s/%s\n', group, name);
fprintf('Size: %d rows x %d columns\n', size(A, 1), size(A, 2));
fprintf('Nonzero entries: %d\n', nnz(A));

%% Display and save the sparsity pattern
figSpy = figure;
spy(A);
title([group '/' name], 'Interpreter', 'none');

exportgraphics(figSpy, imagePath, 'Resolution', 300);
fprintf('Sparsity plot saved to: %s\n', imagePath);

%% Display the matrix as a table
fig = uifigure('Name', [group '/' name]);

uitable(fig, ...
    'Data', full(A), ...
    'Position', [20 20 700 450]);

fig.Position = [100 100 740 490];
