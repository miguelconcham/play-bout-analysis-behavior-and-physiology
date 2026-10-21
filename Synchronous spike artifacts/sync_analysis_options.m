function options = sync_analysis_options(options)
%SYNC_ANALYSIS_OPTIONS  Shared settings for synchronous-spike checks.
%
% Both the waveform pipeline and the raw-LFP power pipeline call this so
% they select the same play-bout coincidences. Pass a struct to override
% any field. Empty fields are filled with the defaults.

if nargin < 1 || isempty(options)
    options = struct();
end

data_root = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data'];

options = fill_option(options, 'data_root', data_root);
options = fill_option(options, 'areas_to_keep', ...
    {'DLPAG', 'LPAG', 'VLPAG', 'DR', 'SupCol'});
options = fill_option(options, 'play_behaviors', ...
    {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', 'Escape', 'CB', 'CD'});
options = fill_option(options, 'hist_range', [-1.5 1]);
options = fill_option(options, 'time_precision', 0.005);
options = fill_option(options, 'min_spikes', 1000);
options = fill_option(options, 'min_sync', 15);
options = fill_option(options, 'n_repeats', 500);
options = fill_option(options, 'min_neuron_fraction', 0.90);
options = fill_option(options, 'max_pairs_per_neuron', inf);
options = fill_option(options, 'max_pairs_per_session', inf);
options = fill_option(options, 'max_pairs', inf);
options = fill_option(options, 'channel_range', 5);
options = fill_option(options, 'n_ap_channels', 384);
options = fill_option(options, 'samples_before', 30);
options = fill_option(options, 'samples_after', 30);
options = fill_option(options, 'n_check_spikes', 40);
options = fill_option(options, 'min_spike_scale', 1);
options = fill_option(options, 'max_first_load_seconds', 600);
options = fill_option(options, 'lfp_sample_rate', 2500);
options = fill_option(options, 'n_lfp_channels', 384);
options = fill_option(options, 'delta_freq', [1 5]);
options = fill_option(options, 'delta_filter_order', 2000);
options = fill_option(options, 'n_phase_bins', 36);
options = fill_option(options, 'phase_guard_seconds', 0.01);
options = fill_option(options, 'raw_folders_file', fullfile(options.data_root, ...
    'NPX data', 'NPX_raw_data_folders.xlsx'));
options = fill_option(options, 'power_definition', ...
    ['square of the Hilbert envelope of unfiltered LFP_PAG, the nearest ' ...
     'sample at each synchronous spike. The null is random LFP samples ' ...
     'inside the same play-bout windows whose 1-5 Hz delta phase falls in ' ...
     'the same 10-degree bin, excluding 10 ms around those spikes']);
end

function options = fill_option(options, name, value)
if ~isfield(options, name) || isempty(options.(name))
    options.(name) = value;
end
end
