function [tgrid, lambda, mean_rate] = kernel_rate(this_session_spikes, t0, t1, dt, sigma)
    % dt = bin size in seconds (e.g. 0.001 or 0.005)
    % sigma = gaussian kernel sd in seconds (e.g. 0.02)
    if nargin < 5 || isempty(sigma), sigma = 0.02; end
    if nargin < 4 || isempty(dt), dt = 0.001; end
    if nargin < 3 || isempty(t0) || isempty(t1)
        if isempty(this_session_spikes)
            t0 = 0; t1 = 1; % fallback
        else
            t0 = min(this_session_spikes);
            t1 = max(this_session_spikes);
        end
    end
    tgrid = t0:dt:t1;
    edges = [tgrid-dt/2, tgrid(end)+dt/2];
    counts = histcounts(this_session_spikes, edges);
    % Gaussian kernel (normalized so integral = 1)
    halfWidth = ceil(5*sigma/dt);
    kx = (-halfWidth:halfWidth)*dt;
    k = exp(-0.5*(kx/sigma).^2);
    k = k / (sum(k)*dt); % normalize so sum(k)*dt = 1 -> convolving returns Hz
    lambda = conv(counts, k, 'same'); % units: spikes per second (Hz)
    mean_rate = mean(lambda); % average of time-varying rate (Hz)
end
