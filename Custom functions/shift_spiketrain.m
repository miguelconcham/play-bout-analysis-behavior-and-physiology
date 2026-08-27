function shifted_spiketrain = shift_spiketrain(spiketrain,starting_time, end_time, time_interval)
% Creates numel(time_interval) shifted spike trains, "circulary" regarding
% starting_time and end_time. The spike train is shifted the aumount of
% time indicated in the time_inteval array each time.

sub_spiketrain = spiketrain(spiketrain<=end_time & spiketrain>=starting_time);

shifted_spiketrain = repmat(sub_spiketrain', numel(time_interval),1) + diag(time_interval)*ones(numel(time_interval), numel(sub_spiketrain));


for interval_index = 1:numel(time_interval)
    
    shifted_spiketrain(interval_index,shifted_spiketrain(interval_index,:)<starting_time) = end_time - starting_time + shifted_spiketrain(interval_index,shifted_spiketrain(interval_index,:)<starting_time);
    shifted_spiketrain(interval_index,shifted_spiketrain(interval_index,:)>end_time) = starting_time - end_time + shifted_spiketrain(interval_index,shifted_spiketrain(interval_index,:)>end_time);
end
end

