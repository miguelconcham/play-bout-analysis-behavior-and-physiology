function Analyze_theta_envelope_delta(result_file)
%ANALYZE_THETA_ENVELOPE_DELTA  Session scan and play-onset test.
%
% Calls Analyze_theta_envelope_session and Analyze_theta_envelope_onset.
% Pass the saved .mat path, or no argument to use the default output of
% Estimate_theta_envelope_delta.

if nargin < 1
    result_file = '';
end
if isempty(result_file)
    Analyze_theta_envelope_session();
    Analyze_theta_envelope_onset();
else
    Analyze_theta_envelope_session(result_file);
    Analyze_theta_envelope_onset(result_file);
end
end
