function token = acg_count_token(method, count_method)
%ACG_COUNT_TOKEN  File-name piece for one autocorrelogram estimate.
%
% Eligible-spike estimates keep _eligible. Overlap pair counts keep the
% original name. FFT and crosscorr insert _fft or _crosscorr, so a new
% count does not replace the pair-count file. A missing count is pairs.

method = canonical_acg_method(method);
if strcmp(method, 'eligible')
    if nargin >= 2 && ~isempty(count_method) ...
            && ~strcmp(canonical_count_method(count_method), 'pairs')
        error('acg_count_token:eligible', ...
            ['Eligible-spike autocorrelograms keep their own count. ' ...
            'acg_count applies to overlap estimates.']);
    end
    token = '_eligible';
    return
end
if nargin < 2 || isempty(count_method)
    count_method = 'pairs';
else
    count_method = canonical_count_method(count_method);
end
switch count_method
    case 'pairs'
        token = '';
    case 'fft'
        token = '_fft';
    case 'crosscorr'
        token = '_crosscorr';
end
end
