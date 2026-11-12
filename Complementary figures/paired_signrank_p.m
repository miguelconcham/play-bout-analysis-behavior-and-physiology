function p = paired_signrank_p(a, b, tail)
%PAIRED_SIGNRANK_P  Two-sided or one-sided signed-rank p-value.
%
% MATLAB signrank is used when it is on the path, which is the test used
% by the locking-period autocorrelogram. Otherwise the normal
% approximation with the tie correction and the continuity correction is
% used. Fewer than two nonzero paired differences returns NaN.

if nargin < 3 || isempty(tail)
    tail = 'both';
end
ok = isfinite(a(:)) & isfinite(b(:));
a = a(ok);
b = b(ok);
if numel(a) < 2 || all(a == b)
    p = NaN;
    return
end
if exist('signrank', 'file') == 2
    try
        p = signrank(a, b, 'tail', tail);
        return
    catch
    end
end
p = wilcoxon_normal_p(a - b, tail);
end


function p = wilcoxon_normal_p(d, tail)
d = d(isfinite(d) & d ~= 0);
n = numel(d);
if n < 2
    p = NaN;
    return
end
ranks = tiedrank(abs(d));
w_pos = sum(ranks(d > 0));
mu = n * (n + 1) / 4;
sigma2 = n * (n + 1) * (2 * n + 1) / 24 - tie_cubic(abs(d)) / 48;
sigma = sqrt(max(sigma2, 0));
if sigma <= 0
    p = NaN;
    return
end
z_of = @(w) continuity_z(w, mu, sigma);
if strcmp(tail, 'right')
    p = normcdf(-z_of(w_pos));
elseif strcmp(tail, 'left')
    p = normcdf(z_of(w_pos));
else
    p = min(1, 2 * normcdf(-abs(z_of(w_pos))));
end
end


function z = continuity_z(w, mu, sigma)
delta = w - mu;
z = (delta - 0.5 * sign(delta)) / sigma;
if delta == 0
    z = 0;
end
end


function s = tie_cubic(v)
sv = sort(v(:));
s = 0;
k = 1;
n = numel(sv);
while k <= n
    j = k;
    while j < n && sv(j + 1) == sv(k)
        j = j + 1;
    end
    t = j - k + 1;
    if t > 1
        s = s + t ^ 3 - t;
    end
    k = j + 1;
end
end
