function obs = speed_peak_pooled_observations(animal, session_key, bin_id, speed, peak)
%SPEED_PEAK_POOLED_OBSERVATIONS  Rows for the pooled peak-frequency models.
%
% SpeedZ is the within-animal z-score of every finite speed bin. The
% returned rows are the comparison set shared by the linear, log-frequency,
% and log-speed models: speed > 0, peak frequency > 0, and a finite
% within-animal z-score. log_peak and log_speed are natural logs.

animal = animal(:);
session_key = session_key(:);
bin_id = double(bin_id(:));
speed = double(speed(:));
peak = double(peak(:));
if numel(peak) ~= numel(speed)
    error('speed_peak_pooled_observations:size', ...
        'peak frequency must have one value per speed.');
end

[speed_z, reference] = within_animal_speed_z(animal, session_key, bin_id, speed);
keep = isfinite(speed) & isfinite(peak) & speed > 0 & peak > 0 & isfinite(speed_z);
obs = struct();
obs.animal = animal(keep);
obs.session_key = session_key(keep);
obs.bin_id = bin_id(keep);
obs.speed = speed(keep);
obs.peak = peak(keep);
obs.log_peak = log(peak(keep));
obs.log_speed = log(speed(keep));
obs.speed_z = speed_z(keep);
obs.n = sum(keep);
obs.n_input = numel(speed);
obs.reference = reference;
end
