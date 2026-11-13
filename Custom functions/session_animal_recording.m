function [animal, recording] = session_animal_recording(session_id)
%SESSION_ANIMAL_RECORDING  Animal and recording from a session folder name.
%
% Session folders look like "B2S2 1110 Single". Animal is the first word.
% Recording is the third word (Single, Dual, or Single2). A name with no
% words keeps the whole string as the animal and leaves the recording blank.

tokens = regexp(strtrim(char(session_id)), '\s+', 'split');
tokens = tokens(~cellfun('isempty', tokens));
if isempty(tokens)
    animal = char(session_id);
    recording = '';
    return
end
animal = tokens{1};
if numel(tokens) >= 3
    recording = tokens{3};
else
    recording = '';
end
end
