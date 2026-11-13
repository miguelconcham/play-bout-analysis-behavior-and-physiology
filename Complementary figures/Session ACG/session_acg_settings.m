function settings = session_acg_settings(options)
%SESSION_ACG_SETTINGS  The fields that decide whether a saved run can be reused.
%
% Window times and neuron ids are not included. Those differ by session.
% bin_edges and smooth_n are stored so a later grand average uses the axis
% the nulls were estimated on.

if ~isfield(options, 'acg_method') || isempty(options.acg_method)
    options.acg_method = 'edge';
end
options.acg_method = char(string(options.acg_method));
if ~any(strcmp(options.acg_method, {'edge', 'crosscorr'}))
    error('session_acg:method', ...
        'acg_method must be ''edge'' or ''crosscorr''.');
end
if ~isfield(options, 'version') || isempty(options.version)
    options.version = 1;
end
required = {'n_shuffles', 'window_s', 'min_spikes', 'tail_percent', ...
    'alpha', 'non_entrained_p', 'oscillation_alpha', 'acg_lag', 'bin_w', ...
    'delta_band', 'smooth_half', 'acg_refractory_s', 'areas_to_keep', ...
    'fooof_max_peaks', 'fooof_fit_range', 'filter_order', 'lfp_sample_rate'};
for i = 1:numel(required)
    if ~isfield(options, required{i}) || isempty(options.(required{i}))
        error('session_acg:settings', 'Missing estimate setting %s.', required{i});
    end
end
if options.n_shuffles < 1 || options.window_s <= 0 || options.min_spikes < 1
    error('session_acg:settings', 'n_shuffles, window_s, and min_spikes must be positive.');
end
if options.tail_percent <= 0 || options.tail_percent >= 50
    error('session_acg:settings', 'tail_percent must be between 0 and 50.');
end
if options.bin_w <= 0 || options.acg_lag <= 0
    error('session_acg:settings', 'acg_lag and bin_w must be positive.');
end
if isfield(options, 'bin_edges') && ~isempty(options.bin_edges)
    bin_edges = double(options.bin_edges(:))';
else
    bin_edges = (-options.acg_lag):options.bin_w:options.acg_lag;
end
settings = struct();
settings.version = double(options.version);
settings.n_shuffles = double(options.n_shuffles);
settings.window_s = double(options.window_s);
settings.min_spikes = double(options.min_spikes);
settings.tail_percent = double(options.tail_percent);
settings.alpha = double(options.alpha);
settings.non_entrained_p = double(options.non_entrained_p);
settings.oscillation_alpha = double(options.oscillation_alpha);
settings.acg_lag = double(options.acg_lag);
settings.bin_w = double(options.bin_w);
settings.delta_band = double(options.delta_band(:))';
settings.smooth_half = double(options.smooth_half);
settings.acg_refractory_s = double(options.acg_refractory_s);
settings.areas_to_keep = cellstr(options.areas_to_keep(:)');
settings.fooof_max_peaks = double(options.fooof_max_peaks);
settings.fooof_fit_range = double(options.fooof_fit_range(:))';
settings.filter_order = double(options.filter_order);
settings.lfp_sample_rate = double(options.lfp_sample_rate);
settings.acg_method = options.acg_method;
settings.bin_edges = bin_edges;
settings.smooth_n = max(1, round(settings.smooth_half / settings.bin_w));
end
