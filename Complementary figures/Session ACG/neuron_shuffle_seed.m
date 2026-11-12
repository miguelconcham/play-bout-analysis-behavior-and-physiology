function seed = neuron_shuffle_seed(session_name, cluster_id)
%NEURON_SHUFFLE_SEED  Stable seed for one neuron's uniform spike times.
%
% The same session and cluster id always redraw the same 500 trains, so a
% resumed run does not need the draws that were already saved.

letters = double(char(string(session_name)));
seed = mod(sum(letters .* (1:numel(letters))) + round(double(cluster_id)) * 997, 2^31 - 2);
seed = seed + 1;
end
