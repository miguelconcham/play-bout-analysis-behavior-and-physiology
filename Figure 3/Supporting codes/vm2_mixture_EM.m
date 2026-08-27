function [mu, kappa, w, LL] = vm2_mixture_EM(theta, maxIter, tol, plotFit)
% VM2_MIXTURE_EM Fit 2-component von Mises mixture using EM
%
% Inputs:
%   theta   - vector of angles in radians (0 to 2*pi)
%   maxIter - maximum iterations (default 100)
%   tol     - convergence tolerance on log-likelihood (default 1e-6)
%   plotFit - if true, plot fitted mixture on circular histogram (default false)
%
% Outputs:
%   mu      - [mu1, mu2] component means
%   kappa   - [kappa1, kappa2] component concentrations
%   w       - weight of first component (second = 1-w)
%   LL      - log-likelihood at each iteration

if nargin < 2, maxIter = 100; end
if nargin < 3, tol = 1e-6; end
if nargin < 4, plotFit = false; end

theta = mod(theta(:), 2*pi);
N = numel(theta);

%% --- Initialization: pick two peaks from KDE ---
grid = linspace(0, 2*pi, 360);
kde = zeros(size(grid));
kappa_kernel = 4; % smoothing
for a = theta'
    kde = kde + exp(kappa_kernel * cos(grid - a));
end
[~, idx] = sort(kde, 'descend');
mu = grid(idx(1:2)); % initial μ1, μ2
kappa = [2, 2];      % initial kappa guesses
w = 0.5;             % equal mixture weight

LL = [];
for iter = 1:maxIter
    %% E-step: responsibilities
    p1 = w * circ_vmpdf(theta, mu(1), kappa(1));
    p2 = (1-w) * circ_vmpdf(theta, mu(2), kappa(2));
    gamma1 = p1 ./ (p1 + p2);
    gamma2 = 1 - gamma1;

    %% M-step: update parameters
    w = mean(gamma1);

    % Component 1
    C1 = sum(gamma1 .* cos(theta));
    S1 = sum(gamma1 .* sin(theta));
    mu(1) = atan2(S1, C1);
    R1 = sqrt(C1^2 + S1^2) / sum(gamma1);
    kappa(1) = circ_kappa(R1);

    % Component 2
    C2 = sum(gamma2 .* cos(theta));
    S2 = sum(gamma2 .* sin(theta));
    mu(2) = atan2(S2, C2);
    R2 = sqrt(C2^2 + S2^2) / sum(gamma2);
    kappa(2) = circ_kappa(R2);

    %% Log-likelihood
    ll = sum(log(w*circ_vmpdf(theta, mu(1), kappa(1)) + ...
                 (1-w)*circ_vmpdf(theta, mu(2), kappa(2))));
    LL = [LL; ll];

    %% Check convergence
    if iter > 1 && abs(LL(end) - LL(end-1)) < tol
        break;
    end
end

%% Optional plot
if plotFit
    figure;
    polarhistogram(theta, 30, 'Normalization','pdf'); hold on;
    theta_grid = linspace(0,2*pi,360);
    pdf_fit = w*circ_vmpdf(theta_grid, mu(1), kappa(1)) + ...
              (1-w)*circ_vmpdf(theta_grid, mu(2), kappa(2));
    polarplot(theta_grid, pdf_fit, 'r','LineWidth',2);
    title('2-component von Mises mixture fit');
end

end

%% --- helper: approximate kappa from R
function k = circ_kappa(R)
% Approximate von Mises kappa from resultant length R
if R < 0.53
    k = 2*R + R^3 + 5*R^5/6;
elseif R < 0.85
    k = -0.4 + 1.39*R + 0.43/(1-R);
else
    k = 1/(R^3 - 4*R^2 + 3*R);
end
k = max(k, 1e-3); % prevent zero
end
