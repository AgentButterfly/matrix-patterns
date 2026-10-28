function A = load_matrix_with_fallback(group, name, dataDir)
%LOAD_MATRIX_WITH_FALLBACK Local MAT/MTX/MTX.GZ, then NIST for BCSSTK/BCSSTM.
% Intended for the real symmetric matrices used in this project.
% Does not use the currently disabled SuiteSparse download endpoint.
%
% Supported automatic archive mapping:
% HB/bcsstk01..13 and HB/bcsstm01..13 -> NIST bcsstruc1
% HB/bcsstk14..18                  -> NIST bcsstruc2
% HB/bcsstk19..25 and HB/bcsstm19..25 -> NIST bcsstruc3
%
% Keep locally cached data: no repeated download if the MAT file exists.

if nargin ~= 3
    error('load_matrix_with_fallback:Arguments', ...
          'Call load_matrix_with_fallback(group,name,dataDir).');
end
if ~isfolder(dataDir), mkdir(dataDir); end
matPath = fullfile(dataDir, [name '.mat']);
mtxPath = fullfile(dataDir, [name '.mtx']);
gzPath  = [mtxPath '.gz'];

if isfile(matPath)
    fprintf('Using existing file: %s\n', matPath);
    loaded = load(matPath, 'Problem');
    if ~isfield(loaded, 'Problem') || ~isfield(loaded.Problem, 'A')
        error('load_matrix_with_fallback:InvalidMAT', ...
            '%s must contain a Problem struct with field A.', matPath);
    end
    A = loaded.Problem.A;
    return
end

if ~isfile(mtxPath)
    if ~isfile(gzPath)
        % Resolve the original archive (only a subset of HB is covered).
        subdir = '';
        matches = regexp(lower(name), '^bcsst[km](\d{2})$', 'tokens', 'once');
        if strcmpi(group,'HB') && ~isempty(matches)
            number = str2double(matches{1});
            if number >= 1 && number <= 13
                subdir = 'bcsstruc1';
            elseif number >= 14 && number <= 18 && ...
                   startsWith(lower(name),'bcsstk')
                subdir = 'bcsstruc2';
            elseif number >= 19 && number <= 25
                subdir = 'bcsstruc3';
            end
        end
        if isempty(subdir)
            error('load_matrix_with_fallback:NoSource', ...
                ['No cached matrix found for %s/%s. SuiteSparse downloads ' ...
                 'are offline. Put %s.mat (containing Problem.A), ' ...
                 '%s.mtx, or %s.mtx.gz in:\n%s\n' ...
                 'This automatic NIST fallback currently supports the ' ...
                 'BCSSTK/BCSSTM Harwell-Boeing series listed in its help.'], ...
                 group,name,name,name,name,dataDir);
        end
        url = sprintf(['https://math.nist.gov/pub/MatrixMarket2/' ...
                       'Harwell-Boeing/%s/%s.mtx.gz'], subdir, name);
        page = sprintf(['https://math.nist.gov/MatrixMarket/data/' ...
                        'Harwell-Boeing/%s/%s.html'],subdir,name);
        fprintf('Downloading %s/%s from NIST: %s\n',group,name,url);
        try
            websave(gzPath, url, weboptions('Timeout',120));
        catch ME
            error('load_matrix_with_fallback:DownloadFailed', ...
                ['Could not download %s/%s from NIST.\n' ...
                 'Open this verified NIST matrix page in a browser:\n%s\n' ...
                 'Download %s.mtx.gz into:\n%s\n' ...
                 'Then rerun the script.\nOriginal error: %s'], ...
                 group,name,page,name,dataDir,ME.message);
        end
    end
    fprintf('Extracting: %s\n', gzPath);
    try
        gunzip(gzPath, dataDir);
    catch ME
        error('load_matrix_with_fallback:UnzipFailed', ...
            'Failed to decompress %s: %s', gzPath, ME.message);
    end
end

fprintf('Reading Matrix Market file: %s\n', mtxPath);
A = read_matrix_market(mtxPath);
Problem = struct('A', A); %#ok<NASGU>
save(matPath, 'Problem', '-v7.3');
fprintf('Saved cached MAT file: %s\n', matPath);
end
