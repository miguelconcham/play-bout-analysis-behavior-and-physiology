function [valsResampled, yResampled] = resample_segments_simple(alignedVals, alignedY, nbins)
% RESAMPLE_SEGMENTS_SIMPLE Resample alignedVals and alignedY into nbins per segment
% following your rules:
%   - If only one 1, fill nbins NaNs before
%   - If two 1s, interpolate between start and first 1
%   - Similarly for 2s

nExp = numel(alignedVals);
valsResampled = cell(1,nExp);
yResampled = cell(1,nExp);

for i = 1:nExp
    y = alignedY{i};
    v = alignedVals{i};
    
    % find indices of 1 and 2 in alignedY
    idx1 = find(y == 1);
    idx2 = find(y == 2);
    
    % --- segment below 1 ---
    if isempty(idx1)
        v0 = nan(1, nbins);
        y0 = linspace(0,1,nbins);
    elseif numel(idx1) == 1
        v0 = nan(1, nbins);
        y0 = linspace(0,1,nbins);
    else
        % two 1s → interpolate values from start to first 1
        v0 = interp1(1:idx1(1), v(1:idx1(1)), linspace(1, idx1(1), nbins), 'linear');
        y0 = linspace(y(1), y(idx1(1)), nbins);
    end
    
    % --- segment 1–2 ---
    if isempty(idx1) || isempty(idx2)
        v1 = nan(1, nbins);
        y1 = linspace(1,2,nbins);
    else
        v1 = interp1(idx1(1):idx2(1), v(idx1(1):idx2(1)), linspace(idx1(1), idx2(1), nbins), 'linear');
        y1 = linspace(y(idx1(1)), y(idx2(1)), nbins);
    end
    
    % --- segment above 2 ---
    if isempty(idx2)
        v2 = nan(1, nbins);
        y2 = linspace(2,3,nbins);
    elseif numel(idx2) == 1
        v2 = nan(1, nbins);
        y2 = linspace(2,3,nbins);
    else
        % two 2s → interpolate from second 2 to last
        v2 = interp1(idx2(2):numel(v), v(idx2(2):end), linspace(idx2(2), numel(v), nbins), 'linear');
        y2 = linspace(y(idx2(2)), y(end), nbins);
    end
    
    % combine segments
    valsResampled{i} = [v0, v1, v2];
    yResampled{i} = [y0, y1, y2];
end
end
