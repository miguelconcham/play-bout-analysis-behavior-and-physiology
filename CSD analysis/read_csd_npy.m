function A = read_csd_npy(filename)
%READ_CSD_NPY  Load a .npy written by write_csd_npy (or compatible NumPy).
%
% Supports little-endian numeric arrays with fortran_order True/False.

fid = fopen(filename, 'rb');
if fid < 0
    error('read_csd_npy:io', 'Could not open %s', filename);
end
cleanup = onCleanup(@() fclose(fid));

magic = fread(fid, 6, 'uint8=>uint8');
if numel(magic) < 6 || ~isequal(magic(:)', uint8([147 78 85 77 80 89]))
    error('read_csd_npy:format', '%s is not a .npy file', filename);
end
major = fread(fid, 1, 'uint8=>uint8');
minor = fread(fid, 1, 'uint8=>uint8'); %#ok<NASGU>
if major == 1
    header_len = fread(fid, 1, 'uint16=>uint16');
elseif major == 2
    header_len = fread(fid, 1, 'uint32=>uint32');
else
    error('read_csd_npy:version', 'Unsupported npy version %d.%d', major, minor);
end
header = char(fread(fid, double(header_len), 'uint8=>char')');

descr = local_parse_descr(header);
fortran = contains(header, '''fortran_order'': True') || contains(header, '"fortran_order": True');
shape = local_parse_shape(header);

switch descr
    case '<f4', precision = 'single=>single';
    case '<f8', precision = 'double=>double';
    case '<i2', precision = 'int16=>int16';
    case '<i4', precision = 'int32=>int32';
    case '<i8', precision = 'int64=>int64';
    otherwise
        error('read_csd_npy:dtype', 'Unsupported descr %s in %s', descr, filename);
end

n = prod(shape);
raw = fread(fid, n, precision);
if numel(raw) ~= n
    error('read_csd_npy:size', 'Expected %d values, got %d in %s', n, numel(raw), filename);
end

if isempty(shape)
    A = raw;
elseif numel(shape) == 1
    A = raw(:);
elseif fortran
    A = reshape(raw, shape);
else
    % C-order on disk: last index varies fastest → reshape then permute
    A = reshape(raw, fliplr(shape));
    A = permute(A, ndims(A):-1:1);
end
end

function descr = local_parse_descr(header)
tok = regexp(header, '''descr'':\s*''([^'']+)''', 'tokens', 'once');
if isempty(tok)
    tok = regexp(header, '"descr":\s*"([^"]+)"', 'tokens', 'once');
end
if isempty(tok)
    error('read_csd_npy:header', 'Could not parse descr from npy header');
end
descr = tok{1};
end

function shape = local_parse_shape(header)
tok = regexp(header, '''shape'':\s*(\([^)]*\))', 'tokens', 'once');
if isempty(tok)
    tok = regexp(header, '"shape":\s*(\([^)]*\))', 'tokens', 'once');
end
if isempty(tok)
    error('read_csd_npy:header', 'Could not parse shape from npy header');
end
inside = tok{1};
inside = strrep(inside, '(', '');
inside = strrep(inside, ')', '');
inside = strtrim(inside);
if isempty(inside)
    shape = [];
    return
end
parts = strsplit(inside, ',');
parts = parts(~cellfun(@(s) isempty(strtrim(s)), parts));
shape = cellfun(@(s) str2double(strtrim(s)), parts);
end
