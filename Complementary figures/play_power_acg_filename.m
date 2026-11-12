function file_name = play_power_acg_filename(animal_name, method)
%PLAY_POWER_ACG_FILENAME  Per-session file for one autocorrelogram method.
%
% Overlap keeps {animal}_play_power_acg.mat. Eligible appends _eligible.
% Passing animal_name '*' is the directory pattern for that method.

method = canonical_acg_method(method);
animal_name = char(animal_name);
if strcmp(method, 'eligible')
    file_name = [animal_name, '_play_power_acg_eligible.mat'];
else
    file_name = [animal_name, '_play_power_acg.mat'];
end
end
