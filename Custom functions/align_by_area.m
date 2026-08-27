function [alignedVals, alignedAreas, alignedY] = align_by_area(allResults, targetAreas)
% Align multiple experiments by placing all available targetAreas into [1,2],
% areas below into [0,1], and areas above into [2,3].

    % Ensure cell array
    if ischar(targetAreas) || isstring(targetAreas)
        targetAreas = {char(targetAreas)};
    end
    
    nExp = numel(allResults);
    alignedVals = {};
    alignedAreas = {};
    alignedY = {};
    
    for iExp = 1:nExp
        for s = 1:numel(allResults(iExp).wrappedVals)
            areas = allResults(iExp).wrappedAreas{s};
            vals  = allResults(iExp).wrappedVals{s};
            y     = allResults(iExp).wrappedY{s};

            if isempty(areas), continue; end

            % --- Step A: detect which target areas are present ---
            isTarget = false(size(areas));
            for a = 1:numel(targetAreas)
                isTarget = isTarget | strcmp(areas, targetAreas{a});
            end
            
            if ~any(isTarget)
                continue; % no relevant areas in this shank
            end

            % --- Step B: compute span of included target areas ---
            yTarget = y(isTarget);
            yMin = min(yTarget);
            yMax = max(yTarget);

            % --- Step C: map into 0–3 normalized axis ---
            yNorm = nan(size(y));
            
            % below block → [0,1]
            belowIdx = y < yMin;
            if any(belowIdx)
                yminAll = min(y(belowIdx));
                ymaxAll = yMin;
                yNorm(belowIdx) = 0 + (y(belowIdx)-yminAll) .* (1-0) ./ (ymaxAll - yminAll);
            end
            
            % target block → [1,2]
            targetIdx = y >= yMin & y <= yMax;
            if any(targetIdx)
                yNorm(targetIdx) = 1 + (y(targetIdx)-yMin) .* (2-1) ./ (yMax - yMin);
            end
            
            % above block → [2,3]
            aboveIdx = y > yMax;
            if any(aboveIdx)
                yminAll = yMax;
                ymaxAll = max(y(aboveIdx));
                yNorm(aboveIdx) = 2 + (y(aboveIdx)-yminAll) .* (3-2) ./ (ymaxAll - yminAll);
            end

            % store
            alignedVals{end+1}  = vals;
            alignedAreas{end+1} = areas;
            alignedY{end+1}     = yNorm;
        end
    end
end
