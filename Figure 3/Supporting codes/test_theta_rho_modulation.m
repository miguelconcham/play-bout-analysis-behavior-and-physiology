function [pval, results] = test_theta_rho_modulation(theta, rho, nbins, nperm, mode)
%TEST_THETA_RHO_MODULATION  Test whether rho depends on theta.
%
%   [pval, results] = test_theta_rho_modulation(theta, rho, nbins, nperm, mode)
%
%   Inputs:
%     theta : vector of angles (radians)
%     rho   : vector of associated values
%     nbins : number of angular bins (default: 12)
%     nperm : number of permutations (default: 1000)
%     mode  : 'modulation' (default) or 'increase'
%
%   Outputs:
%     pval     : permutation test p-value
%     results  : struct with fields
%                   .bin_centers
%                   .bin_means
%                   .observed_stat   (depends on mode)
%                   .null_stats      (permutation null distribution)
%                   .theta_max       (bin center with max rho)
%                   .rho_max         (max mean rho)
%                 if mode='modulation':
%                   .theta_min       (bin center with min rho)
%                   .rho_min         (min mean rho)

    if nargin < 3 || isempty(nbins)
        nbins = 12;
    end
    if nargin < 4 || isempty(nperm)
        nperm = 1000;
    end
    if nargin < 5 || isempty(mode)
        mode = 'modulation'; % default
    end
    
    theta = mod(theta, 2*pi); % wrap to [0, 2pi)
    
    % Bin edges and centers
    edges = linspace(0, 2*pi, nbins+1);
    bin_centers = edges(1:end-1) + diff(edges)/2;
    
    % Compute mean rho in each bin
    bin_means = nan(1, nbins);
    for b = 1:nbins
        idx = theta >= edges(b) & theta < edges(b+1);
        bin_means(b) = mean(rho(idx), 'omitnan');
    end
    
    % Observed statistics
    [rho_max, idx_max] = max(bin_means);
    results.theta_max = bin_centers(idx_max);
    results.rho_max   = rho_max;
    
    switch lower(mode)
        case 'modulation'
            [rho_min, idx_min] = min(bin_means);
            observed_stat = rho_max - rho_min;
            results.theta_min = bin_centers(idx_min);
            results.rho_min   = rho_min;
            
        case 'increase'
            observed_stat = rho_max;
            
        otherwise
            error('Unknown mode: %s. Use ''modulation'' or ''increase''.', mode);
    end
    
    % Permutation null distribution
    null_stats = nan(1, nperm);
    for ip = 1:nperm
        rho_perm = rho(randperm(numel(rho)));
        tmp_means = nan(1, nbins);
        for b = 1:nbins
            idx = theta >= edges(b) & theta < edges(b+1);
            tmp_means(b) = mean(rho_perm(idx), 'omitnan');
        end
        
        switch lower(mode)
            case 'modulation'
                null_stats(ip) = max(tmp_means) - min(tmp_means);
            case 'increase'
                null_stats(ip) = max(tmp_means);
        end
    end
    
    % p-value
    pval = mean(null_stats >= observed_stat);
    
    % Package results
    results.bin_centers   = bin_centers;
    results.bin_means     = bin_means;
    results.observed_stat = observed_stat;
    results.null_stats    = null_stats;
end
