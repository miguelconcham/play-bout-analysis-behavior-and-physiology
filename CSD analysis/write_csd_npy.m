function write_csd_npy(filename, A)
%WRITE_CSD_NPY  Write a numeric array as .npy readable by NumPy.
%
% Matches the common npy-matlab / writeNPY convention: Fortran order, so
% np.load(filename) has the same shape as size(A) in MATLAB.

if nargin < 2
    error('write_csd_npy:args', 'Usage: write_csd_npy(filename, A)');
end

if ~isnumeric(A) || isempty(A)
    error('write_csd_npy:type', 'A must be a non-empty numeric array.');
end

if isa(A, 'single')
    dtype = 'f4';
elseif isa(A, 'double')
    dtype = 'f8';
elseif isa(A, 'int16')
    dtype = 'i2';
elseif isa(A, 'int32')
    dtype = 'i4';
elseif isa(A, 'int64')
    dtype = 'i8';
else
    A = double(A);
    dtype = 'f8';
end

shape = size(A);
if iscolumn(A)
    shape_np = numel(A);
elseif isrow(A)
    shape_np = [1, numel(A)];
else
    shape_np = shape;
end

% Build shape text as char (avoid string arrays — fwrite rejects them)
parts = cell(1, numel(shape_np));
for k = 1:numel(shape_np)
    parts{k} = num2str(shape_np(k));
end
if numel(shape_np) == 1
    shape_str = ['(' parts{1} ',)'];
else
    shape_str = ['(' strjoin(parts, ', ') ')'];
end
shape_str = char(shape_str);
dtype = char(dtype);

header_dict = sprintf( ...
    '{''descr'': ''<%s'', ''fortran_order'': True, ''shape'': %s, }', ...
    dtype, shape_str);
header_dict = char(header_dict);

% npy v1.0: magic(6) + ver(2) + header_len(uint16) = 10-byte prefix.
% Header + prefix must be aligned; use 64 bytes (stricter than required 16).
prefix_len = 10;
pad_len = 64 - mod(prefix_len + numel(header_dict) + 1, 64);
if pad_len == 64
    pad_len = 0;
end
header = [header_dict, repmat(char(' '), 1, pad_len), char(10)];
header = char(header(:).');  % 1 x N char row
header_len = uint16(numel(header));

fid = fopen(filename, 'wb');
if fid < 0
    error('write_csd_npy:io', 'Could not open %s for writing', filename);
end
fwrite(fid, uint8([147 78 85 77 80 89]), 'uint8'); % \x93NUMPY
fwrite(fid, uint8(1), 'uint8');  % major
fwrite(fid, uint8(0), 'uint8');  % minor
fwrite(fid, header_len, 'uint16');
fwrite(fid, uint8(header), 'uint8');  % safe binary write (not 'char'/string)
fwrite(fid, A(:), class(A));
fclose(fid);
end
