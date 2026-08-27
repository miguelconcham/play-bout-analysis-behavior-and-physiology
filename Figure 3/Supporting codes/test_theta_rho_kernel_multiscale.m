function [pval, results] = test_theta_rho_kernel_multiscale(theta, rho, kappa_list, nperm, mode, ngr)
%TEST_THETA_RHO_KERNEL_MULTISCALE  Multiscale circular-kernel test.
%
% [pval, results] = test_theta_rho_kernel_multiscale(theta, rho, kappa_list, nperm, mode, ngr)
%
% Inputs:
%  - theta: angles (rad), can be in [-pi,pi] or [0,2pi)
%  - rho:   corresponding values (same length)
%  - kappa_list: vector of von Mises concentration params (e.g. [1 5 20 80])
%       larger kappa -> narrower kernel
%  - nperm: number of permutations (default 2000)
%  - mode: 'increase' (default) or 'modulation' (two-sided)
%  - ngr: number of evaluation grid points around circle (default 360)
%
% Outputs:
%  - pval: permutation p-value (one-sided for 'increase' or two-sided style if 'modulation')
%  - results: struct with fields
%       .theta_grid (ngr points)
%       .kappa_list
%       .smoothed_obs{ik} (smoothed rho over theta_grid for each kappa)
%       .observed_stat (global max across kappa/grid or max-min if modulation)
%       .peak_theta (angle where max observed occurs, in same range as input)
%       .peak_kappa (kappa producing that peak)
%       .null_stats (nperm null global stats)
%
% Example:
%  [p, res] = test_theta_rho_kernel_multiscale(theta, rho, [2 8 32], 1000, 'increase');
%

if nargin < 3 || isempty(kappa_list), kappa_list = [2 8 32]; end
if nargin < 4 || isempty(nperm), nperm = 2000; end
if nargin < 5 || isempty(mode), mode = 'increase'; end
if nargin < 6 || isempty(ngr), ngr = 360; end

% Keep input theta range to return results in same range
theta_in = theta;
wrap_to_pi_flag = any(theta < 0); % if input had negative values assume [-pi,pi] convention

% normalize angles to [0,2pi)
theta = mod(theta, 2*pi);
N = numel(theta);

% evaluation grid
theta_grid = linspace(0, 2*pi, ngr+1);
theta_grid(end) = []; % remove duplicate 2pi
K = numel(kappa_list);

% pre-allocate
smoothed_obs = cell(1,K);
obs_stats = zeros(1,K);

% von Mises kernel smoother at each kappa
for ik = 1:K
    kappa = kappa_list(ik);
    % compute smoothed on grid
    S = zeros(1,ngr);
    for ig = 1:ngr
        delta = theta - theta_grid(ig);
        w = exp(kappa * cos(delta));   % von Mises kernel (unnormalized)
        S(ig) = sum(w .* rho) / sum(w);
    end
    smoothed_obs{ik} = S;
    switch lower(mode)
        case 'increase'
            obs_stats(ik) = max(S);           % max over grid
        case 'modulation'
            obs_stats(ik) = max(S) - min(S);  % peak-to-trough
        otherwise
            error('Unknown mode. Use ''increase'' or ''modulation''.');
    end
end

% observed global stat (max across kappas)
[observed_stat, idxK] = max(obs_stats);
% find grid index and angle of peak for 'increase'. For 'modulation', report theta of max.
if strcmpi(mode,'increase')
    Sbest = smoothed_obs{idxK};
    [~, ig_peak] = max(Sbest);
    peak_theta = theta_grid(ig_peak);
else
    Sbest = smoothed_obs{idxK};
    [~, ig_peak] = max(Sbest);
    peak_theta = theta_grid(ig_peak);
end

% Permutation null distribution: compute global stat for each perm
null_stats = nan(1, nperm);
rng('shuffle'); % randomize seed for permutations
for ip = 1:nperm
    rho_perm = rho(randperm(N));
    tmp_stats = zeros(1,K);
    for ik = 1:K
        kappa = kappa_list(ik);
        Sperm = zeros(1,ngr);
        for ig = 1:ngr
            delta = theta - theta_grid(ig);
            w = exp(kappa * cos(delta));
            Sperm(ig) = sum(w .* rho_perm) / sum(w);
        end
        if strcmpi(mode,'increase')
            tmp_stats(ik) = max(Sperm);
        else
            tmp_stats(ik) = max(Sperm) - min(Sperm);
        end
    end
    null_stats(ip) = max(tmp_stats); % global max across kappas
end

% p-value (one-sided tail: fraction null >= observed)
pval = mean(null_stats >= observed_stat);

% pack results, convert peak angle back to input convention
if wrap_to_pi_flag
    peak_theta_out = wrapToPi(peak_theta);
else
    peak_theta_out = peak_theta;
end

results.theta_grid = theta_grid;
results.kappa_list = kappa_list;
results.smoothed_obs = smoothed_obs;
results.observed_stat = observed_stat;
results.peak_theta = peak_theta_out;
results.peak_kappa = kappa_list(idxK);
results.null_stats = null_stats;
results.mode = mode;

end
