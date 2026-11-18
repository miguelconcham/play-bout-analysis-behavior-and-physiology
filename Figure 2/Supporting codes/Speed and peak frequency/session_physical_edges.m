function [physical_edges, found] = session_physical_edges(psth_structure, session_rows, z_name, physical_name, z_edges)
%SESSION_PHYSICAL_EDGES  Convert each session's z edge with its own scale.
%
% A within-session z of 1 is a different number of centimetres in every
% session. The physical column must be paired with that session's z-score.

session_rows = session_rows(:);
z_edges = z_edges(:);
n_sess = numel(session_rows);
physical_edges = nan(n_sess, 1);
found = false;
if numel(z_edges) ~= n_sess
    error('session_physical_edges:Size', ...
        'Expected one z edge for each of %d sessions.', n_sess);
end
for s = 1:n_sess
    model_data = psth_structure(session_rows(s)).model_data;
    [z_samples, physical_samples] = aligned_z_and_physical(model_data, z_name, physical_name);
    if isempty(physical_samples)
        continue
    end
    found = true;
    physical_edges(s) = z_to_physical(z_edges(s), z_samples, physical_samples);
end
end
