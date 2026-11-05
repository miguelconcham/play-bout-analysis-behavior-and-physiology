function options = irasa_options(options)
%IRASA_OPTIONS  Settings for irregular-resampling auto-spectral analysis.
%
% Wen & Liu (2016), Brain Topography 29:13-26, section 2.4. Pass a struct
% to override any field. Empty fields keep the default.
%
% h_values are the resampling factors h > 1. Each h is paired with 1/h.
% The published set is 1.1, 1.15, ..., 1.9.
% n_segments and segment_fraction are the ten windows, each 90% of the
% epoch, in section 2.4 step 1.
% detrend removes a least-squares line from each window before the FFT.
% antialias low-passes before the cubic-spline downsampling in step 3.
% frequency_range is [fmin fmax] in Hz. Leave it empty to keep every
% positive frequency up to the bandwidth that survives the largest h,
% which is sampling_rate / (2 * max(h_values)).

if nargin < 1 || isempty(options)
    options = struct();
end

options = fill_option(options, 'h_values', 1.1:0.05:1.9);
options = fill_option(options, 'n_segments', 10);
options = fill_option(options, 'segment_fraction', 0.9);
options = fill_option(options, 'detrend', true);
options = fill_option(options, 'antialias', true);
options = fill_option(options, 'frequency_range', []);

validate_options(options);
end

function options = fill_option(options, name, value)
if ~isfield(options, name) || isempty(options.(name))
    options.(name) = value;
end
end

function validate_options(options)
h_values = options.h_values(:);
if ~isnumeric(h_values) || ~all(isfinite(h_values)) || any(h_values <= 1)
    error('irasa_options:h', ...
        'h_values must be finite resampling factors greater than 1.');
end
if ~isscalar(options.n_segments) || options.n_segments < 1 || ...
        options.n_segments ~= round(options.n_segments)
    error('irasa_options:segments', ...
        'n_segments must be a positive integer.');
end
if ~isscalar(options.segment_fraction) || ~isfinite(options.segment_fraction) || ...
        options.segment_fraction <= 0 || options.segment_fraction > 1
    error('irasa_options:fraction', ...
        'segment_fraction must lie in (0, 1].');
end
if options.n_segments > 1 && options.segment_fraction == 1
    error('irasa_options:overlap', ...
        'Overlapping segments need segment_fraction below 1.');
end
if ~(islogical(options.detrend) || isnumeric(options.detrend)) || ...
        ~isscalar(options.detrend)
    error('irasa_options:detrend', 'detrend must be a scalar logical.');
end
if ~(islogical(options.antialias) || isnumeric(options.antialias)) || ...
        ~isscalar(options.antialias)
    error('irasa_options:antialias', 'antialias must be a scalar logical.');
end
frequency_range = options.frequency_range;
if ~isempty(frequency_range)
    if ~isnumeric(frequency_range) || numel(frequency_range) ~= 2 || ...
            ~all(isfinite(frequency_range)) || frequency_range(1) < 0 || ...
            frequency_range(2) <= frequency_range(1)
        error('irasa_options:frequency', ...
            'frequency_range must be empty or [fmin fmax] with fmax > fmin.');
    end
end
end
