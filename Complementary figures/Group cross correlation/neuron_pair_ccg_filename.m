function name = neuron_pair_ccg_filename(animal_code)
%NEURON_PAIR_CCG_FILENAME  Per-session neuron-pair file.

name = [ccg_plain_text(animal_code), '_neuron_pair_ccg.mat'];
end
