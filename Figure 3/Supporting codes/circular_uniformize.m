function uniform_phases = circular_uniformize(current_angles)
% CIRCULAR_UNIFORMIZE  Transform circular data into a uniform phase distribution.
%
%   uniform_phases = circular_uniformize(current_angles)
%
%   This function maps an arbitrary angular distribution (e.g. LFP phases)
%   onto a uniform distribution on the circle [-pi, pi]. It does so by:
%     1. Wrapping all input angles into [-pi, pi].
%     2. Sorting angles and assigning empirical cumulative probabilities.
%     3. Extending the domain by ±2π to handle circular wrap-around.
%     4. Interpolating the CDF at the original angles.
%     5. Scaling the uniform values to [-pi, pi].
%
%   INPUT:
%       current_angles - Vector of angles (in radians). Can be any length.
%                        Values will be wrapped into [-pi, pi].
%
%   OUTPUT:
%       uniform_phases - Vector of the same size as input, with approximately
%                        uniform distribution on [-pi, pi].
%
%   NOTES:
%   - This is the circular equivalent of applying the probability integral
%     transform (CDF transform) to linear variables.
%   - Useful for spike–phase coupling analysis where you want to remove
%     bias due to non-uniform LFP phase distributions.
%
%   See also: ANGLE, MOD, INTERP1
    % Ensure angles are in [-pi, pi]
    current_angles = mod(current_angles+pi, 2*pi) - pi;

    % Sort angles
    [sorted_angles, idx] = sort(current_angles);

    % Circular empirical CDF (uniform spacing around circle)
    f = ((1:numel(sorted_angles)) - 0.5) / numel(sorted_angles);

    % Map back to original angles using interpolation (circular wrap)
    % Extend domain by wrapping endpoints
    x_ext = [sorted_angles-2*pi, sorted_angles, sorted_angles+2*pi];
    f_ext = [f-1, f, f+1];   % wrap CDF accordingly

    [x_ext, ia] = unique(x_ext, 'stable');
    f_ext = f_ext(ia);
    % Interpolate each original angle
    u = interp1(x_ext, f_ext, current_angles, 'linear');

    % Now u ~ Uniform[0,1], map to [-pi, pi]
    uniform_phases = 2*pi*u - pi;
end