function recordings = GENERATE_THETA_ENVELOPE_DELTA(current_dir, options)
%GENERATE_THETA_ENVELOPE_DELTA  Theta-envelope delta estimates for one recording.
%
% One struct per mid-LPAG channel. Session fields hold the envelope spectrum,
% its coherence with the LFP, and surrogate z-scores on the 0.5 Hz and 0.2 Hz
% grids. Onset fields hold play-aligned epochs and the onset-shift null.

if nargin < 2 || isempty(options)
    options = struct();
end
options = theta_envelope_delta_options(options);
options.theta_coefficients = theta_envelope_filter_coefficients(options);
options.grid_test = theta_envelope_frequency_grid( ...
    options.modulation_band, options.test_step_hz, options.bandwidth_hz, ...
    options.confirmatory_band);
options.grid_fine = theta_envelope_frequency_grid( ...
    options.modulation_band, options.fine_step_hz, options.bandwidth_hz, ...
    options.confirmatory_band);

loaded = theta_envelope_load_recording(current_dir, options);
recordings = [];
if isempty(loaded.selected)
    return
end
for channel_n = 1:numel(loaded.selected)
    selected = loaded.selected(channel_n);
    fprintf('%s  probe %g  channel %g\n', loaded.animal_code, ...
        selected.probe_number, selected.channel);
    one = theta_envelope_one_channel( ...
        loaded.lfp(channel_n, :), selected, loaded, options);
    if isempty(one)
        continue
    end
    if isempty(recordings)
        recordings = one;
    else
        recordings(end + 1) = one; %#ok<AGROW>
    end
end
end
