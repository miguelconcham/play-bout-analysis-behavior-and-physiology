function settings = group_ccg_settings(options)
%GROUP_CCG_SETTINGS  Lag axis, areas, and the minimum reference count.
%
% The estimate bin is 1 ms. Edges run from -acg_lag to +acg_lag, and the
% step is acg_lag divided by the number of steps so the ends land on the
% requested lag. Zero is an edge, so the two bins that meet at zero split
% simultaneous spikes into the non-negative bin. areas_to_keep are the
% areas whose spikes are pooled. SupCol is accepted as SC by
% canonical_kept_area.

if nargin < 1 || isempty(options)
    options = struct();
end
settings = struct();
settings.version = 2;
settings.acg_lag = option_or(options, 'acg_lag', 1);
settings.acg_bin = option_or(options, 'acg_bin', 0.001);
settings.areas_to_keep = option_or(options, 'areas_to_keep', ...
    {'SC', 'LPAG', 'VLPAG', 'DLPAG', 'DR'});
if isfield(options, 'session_type')
    settings.session_type = ccg_plain_text(options.session_type);
else
    settings.session_type = 'Partners session';
end
settings.min_reference_spikes = option_or(options, 'min_reference_spikes', 50);
settings.alpha = option_or(options, 'alpha', 0.01);
settings.non_entrained_p = option_or(options, 'non_entrained_p', 0.1);
if ~(settings.acg_lag > 0) || ~(settings.acg_bin > 0)
    error('group_ccg:settings', 'acg_lag and acg_bin must be positive.');
end
if settings.acg_bin >= settings.acg_lag
    error('group_ccg:settings', 'acg_bin must be shorter than acg_lag.');
end
if settings.min_reference_spikes < 1
    error('group_ccg:settings', 'min_reference_spikes must be at least 1.');
end
if settings.alpha <= 0 || settings.non_entrained_p <= settings.alpha
    error('group_ccg:settings', ...
        'non_entrained_p must be greater than alpha, and alpha must be positive.');
end
settings.areas_to_keep = cellstr(settings.areas_to_keep(:))';
settings.session_type = ccg_plain_text(settings.session_type);
n_side = round(settings.acg_lag / settings.acg_bin);
if abs(n_side * settings.acg_bin - settings.acg_lag) > 1e-6 * settings.acg_lag
    error('group_ccg:settings', 'acg_lag must be an integer multiple of acg_bin.');
end
settings.acg_bin = settings.acg_lag / n_side;
settings.bin_edges = (-n_side:n_side) * settings.acg_bin;
settings.centers = settings.bin_edges(1:end - 1) + settings.acg_bin / 2;
end


function value = option_or(options, name, default_value)
if isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = default_value;
end
end
