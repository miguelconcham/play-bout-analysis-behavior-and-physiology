function file_name = probe_periodic_spectra_filename(animal_name)
% File name shared by the channel-wise partner-session spectra.
%
% A Neuropixels folder named "B1D1 1007 Dual" is saved as
% play_session_probe_periodic_spectra_B1D1_1007_Dual.mat.
% Spaces in the folder name become single underscores. The three parts are
% the animal, the date, and Single, Dual, or Single2.
%
% file_name = probe_periodic_spectra_filename()
%   Dir pattern for every recording saved under this stem.
% file_name = probe_periodic_spectra_filename(animal_name)
%   File name for one recording folder.

stem = 'play_session_probe_periodic_spectra_';
if nargin < 1 || isempty(animal_name)
    file_name = {[stem, '*']};
    return
end
session_name = strtrim(char(animal_name));
session_name = regexprep(session_name, '\s+', '_');
file_name = [stem, session_name, '.mat'];
end
