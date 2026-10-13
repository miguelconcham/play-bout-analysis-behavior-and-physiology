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
    shape_np = [numel(A)];
elseif isrow(A)
    shape_np = [1, numel(A)];
else
    shape_np = shape;
end

shape_str = ['(' strjoin(arrayfun(@num2str, shape_np, 'UniformOutput', false), ', ')];
if numel(shape_np) == 1
    shape_str = [shape_str, ',)'];
else
    shape_str = [shape_str, ')'];
end

header_dict = sprintf( ...
    "{'descr': '<%s', 'fortran_order': True, 'shape': %s, }", ...
    dtype, shape_str);

magic = [147 78 85 77 80 89]; % \x93NUMPY
major = uint8(1);
minor = uint8(0);
header_len_pos = 10;
pad_len = 64 - mod(header_len_pos + 2 + numel(header_dict) + 1, 64);
if pad_len == 64
    pad_len = 0;
end
header = [header_dict, repmat(' ', 1, pad_len), sprintf('\n')];
header_len = uint16(numel(header));

fid = fopen(filename, 'wb');
if fid < 0
    error('write_csd_npy:io', 'Could not open %s for writing', filename);
end
fwrite(fid, magic, 'uint8');
fwrite(fid, major, 'uint8');
fwrite(fid, minor, 'uint8');
fwrite(fid, header_len, 'uint16');
fwrite(fid, header, 'char');
fwrite(fid, A(:), class(A));
fclose(fid);
end
