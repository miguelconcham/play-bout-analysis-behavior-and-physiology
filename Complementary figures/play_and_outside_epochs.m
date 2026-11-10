function epochs = play_and_outside_epochs(session_windows, play_bouts, epoch_duration, time_limits)
%PLAY_AND_OUTSIDE_EPOCHS  Fixed epochs inside play and outside play.
%
% session_windows and play_bouts are N x 2, [start, end] in seconds.
% time_limits is [first, last] and is optional. Play is first restricted
% to the session windows, and both are clipped to time_limits when those
% are given, so the epochs stay inside the LFP. Epochs are
% epoch_duration long and never cross a play boundary.

if nargin < 4
    time_limits = [];
end
if ~isempty(time_limits)
    session_windows = intersect_time_intervals(session_windows, time_limits);
    play_bouts = intersect_time_intervals(play_bouts, time_limits);
end
play_bouts = intersect_time_intervals(play_bouts, session_windows);
epochs.play_bouts = play_bouts;
epochs.outside_intervals = outside_play_intervals(session_windows, play_bouts);
epochs.play = fixed_duration_epochs(play_bouts, epoch_duration);
epochs.outside = fixed_duration_epochs(epochs.outside_intervals, epoch_duration);
end
