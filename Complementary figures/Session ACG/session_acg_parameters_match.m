function match = session_acg_parameters_match(left, right)
%SESSION_ACG_PARAMETERS_MATCH  True when two estimate settings describe one run.

match = isequaln(session_acg_settings(left), session_acg_settings(right));
end
