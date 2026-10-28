function options = playbout_erp_defaults()
%PLAYBOUT_ERP_DEFAULTS  Shared settings for the play-bout raw-LFP ERP.
%
% The epoch is the full unfiltered mid-LPAG waveform from -2 s to +2 s
% around play-bout onset. Saturated trials are those that reach the int16
% rails of LFP_PAG. The analysis baseline is the pre-onset half of that
% epoch, [−2, 0).

options.window_s = [-2 2];
options.sampling_rate = 2500;
options.saturation_limit = 32767;
options.baseline_window = [-2 0];
options.electrode_index = 1;
options.lme_step_samples = 1;
options.alpha = 0.05;
options.formula = 'Voltage ~ 1 + (1|Animal)';
options.play_behaviors = {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', ...
    'Escape', 'CB', 'CD'};
options.area = 'PAG';
options.lfp_units = 'LFP_PAG sample values';
end
