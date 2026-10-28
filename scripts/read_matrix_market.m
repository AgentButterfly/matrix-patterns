function A = read_matrix_market(filePath)
%READ_MATRIX_MARKET Read a real coordinate Matrix Market .mtx file.
% Supports general and symmetric real/integer coordinate matrices.
% The full symmetric matrix is reconstructed when the file stores only
% one triangular part. Intended for the first HB pattern-analysis tests.

fid = fopen(filePath, 'rt');
if fid < 0
    error('read_matrix_market:OpenFailed', 'Cannot open: %s', filePath);
end
cleanup = onCleanup(@() fclose(fid));

header = strtrim(fgetl(fid));
tokens = strsplit(lower(header));
if numel(tokens) ~= 5 || ...
        ~strcmp(tokens{1}, '%%matrixmarket') || ...
        ~strcmp(tokens{2}, 'matrix') || ...
        ~strcmp(tokens{3}, 'coordinate') || ...
        ~ismember(tokens{4}, {'real','integer'}) || ...
        ~ismember(tokens{5}, {'general','symmetric'})
    error('read_matrix_market:UnsupportedFormat', ...
        'Expected real/integer coordinate, general/symmetric Matrix Market: %s', ...
        header);
end
isSymmetric = strcmp(tokens{5}, 'symmetric');

% Skip comment lines to read the size declaration.
line = fgetl(fid);
while ischar(line) && (isempty(strtrim(line)) || ...
                       startsWith(strtrim(line), '%'))
    line = fgetl(fid);
end
if ~ischar(line)
    error('read_matrix_market:InvalidHeader', 'Missing matrix dimensions.');
end
sizes = sscanf(line, '%d %d %d');
if numel(sizes) ~= 3
    error('read_matrix_market:InvalidHeader', 'Invalid matrix dimensions.');
end
m = sizes(1); n = sizes(2); nz = sizes(3);
if isSymmetric && m ~= n
    error('read_matrix_market:InvalidHeader', ...
        'Symmetric Matrix Market input must be square.');
end

values = textscan(fid, '%f %f %f', nz, 'CollectOutput', true);
triples = values{1};
if size(triples,1) ~= nz || size(triples,2) ~= 3
    error('read_matrix_market:InvalidEntries', ...
        'Expected %d coordinate entries.', nz);
end
i = triples(:,1); j = triples(:,2); v = triples(:,3);
if any(i < 1 | i > m | j < 1 | j > n | ...
       i ~= floor(i) | j ~= floor(j))
    error('read_matrix_market:InvalidEntries', 'Invalid coordinate indices.');
end
A = sparse(i, j, v, m, n);
if isSymmetric
    % For Matrix Market symmetric format, off-diagonal entries are stored
    % only once. Add their mirrored entries, but count the diagonal once.
    A = A + A.' - spdiags(diag(A), 0, n, n);
end
end
