function ppc = compute_ppc_fast(phases)
% COMPUTE_PPC_FAST Efficient PPC estimate for large spike counts
%
%   ppc = compute_ppc_fast(phases)
%
%   Input:
%       phases : vector of spike phases in radians
%   Output:
%       ppc    : scalar PPC value

    phases = phases(:);           % ensure column vector
    N = numel(phases);
    
    if N < 2
        warning('Need at least 2 spikes to compute PPC.');
        ppc = NaN;
        return;
    end
    
    % Compute complex phase vectors
    z = exp(1i * phases);
    
    % Efficient PPC calculation
    R = sum(z);
    ppc = (abs(R)^2 - N) / (N * (N - 1));
end
