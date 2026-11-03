function choice = speed_peak_pooled_formula(animal, session_key)
%SPEED_PEAK_POOLED_FORMULA  Repeated-measure formula for the pooled models.
%
% The intended model is outcome ~ predictor + (1|Animal) + (1|Session).
% Session is nested in animal. With one session per animal the two
% intercepts are the same term, so the formula keeps the animal intercept.
% With one animal the animal intercept is the fixed intercept and sessions
% carry the repeated measures.

animal = animal(:);
session_key = session_key(:);
if ~iscell(animal) || ~iscell(session_key) || numel(animal) ~= numel(session_key)
    error('speed_peak_pooled_formula:size', ...
        'animal and session_key must be cell arrays of the same length.');
end

n_animals = numel(unique(animal));
n_sessions = numel(unique(session_key));
choice = struct( ...
    'n_animals', n_animals, ...
    'n_sessions', n_sessions, ...
    'random_terms', '', ...
    'note', '');
if isempty(animal)
    choice.note = 'No rows.';
    return
end
if n_animals >= 2 && n_sessions > n_animals
    choice.random_terms = '(1|Animal) + (1|Session)';
    choice.note = ['Repeated speed bins within session, and repeated ' ...
        'sessions within animal.'];
elseif n_animals >= 2
    choice.random_terms = '(1|Animal)';
    choice.note = ['Each animal contributes one session, so the session ' ...
        'intercept is the animal intercept.'];
elseif n_sessions >= 2
    choice.random_terms = '(1|Session)';
    choice.note = ['One animal. The repeated measures are its sessions, ' ...
        'and the animal intercept is the fixed intercept.'];
else
    choice.note = ['One animal and one session. The fit is ordinary least ' ...
        'squares on the repeated speed bins.'];
end
end
