%% Analyze_theta_envelope_delta
% Section script. Run it from the top, or run one section at a time.
%
% %% Session spectrum and coherence
%     opens the envelope-power figure and the envelope-LFP coherence figure
% %% Play onset
%     opens, for the envelope and for spectrogram power, and for the
%     implanted animal and the partner: one figure with a panel per
%     behavior, and one figure that stacks those behaviors above their mean.
%
% The two sections below only call the other section scripts. Open
% Analyze_theta_envelope_session.m or Analyze_theta_envelope_onset.m to run
% a single figure.

%% Paths
this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Figure 2\' ...
        'Supporting codes v2\Theta envelope\Analyze_theta_envelope_delta.m'];
end
repo_root = fileparts(fileparts(fileparts(fileparts(this_file))));
run(fullfile(repo_root, 'add_repo_paths.m'));

%% Session spectrum and coherence
Analyze_theta_envelope_session

%% Play onset
Analyze_theta_envelope_onset
