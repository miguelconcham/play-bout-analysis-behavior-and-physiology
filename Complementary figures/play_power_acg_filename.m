function file_name = play_power_acg_filename(animal_name, method, count_method)
%PLAY_POWER_ACG_FILENAME  Per-session file for one autocorrelogram method.
%
% Overlap pair counts keep {animal}_play_power_acg.mat. Eligible appends
% _eligible. FFT and crosscorr append _fft or _crosscorr. Passing
% animal_name '*' is the directory pattern for that method and count.

if nargin < 3
    count_method = [];
end
token = acg_count_token(method, count_method);
animal_name = char(animal_name);
file_name = [animal_name, '_play_power_acg', token, '.mat'];
end
