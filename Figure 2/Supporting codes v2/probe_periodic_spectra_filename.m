function file_name = probe_periodic_spectra_filename(animal_name)
% File name shared by the channel-wise partner-session spectra.
%
% file_name = probe_periodic_spectra_filename()
%   Cell of dir patterns. The current stem is
%   play_session_probe_periodic_spectra_. The earlier stem
%   play_session_probe_periodic_aperiodic_spectra_ is included so files
%   saved under either name are found without editing the scripts.
% file_name = probe_periodic_spectra_filename(animal_name)
%   Name used when saving one recording. Separators in animal_name become
%   underscores.

current_stem = 'play_session_probe_periodic_spectra_';
if nargin < 1 || isempty(animal_name)
    file_name = { ...
        [current_stem, '*.mat'], ...
        'play_session_probe_periodic_aperiodic_spectra_*.mat'};
    return
end
safe_name = regexprep(char(animal_name), '[^\w\-]+', '_');
file_name = [current_stem, safe_name, '.mat'];
end
