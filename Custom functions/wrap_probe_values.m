function [wrappedVals, wrappedAreas, wrappedY, borders] = ...
    wrap_probe_values(ARRAY1, CELLARRAY2, ARRAY3, nbins, probeType)
% ARRAY1: [nchann x 2] (x,y)
% CELLARRAY2: {nchann x 1} area names
% ARRAY3: [nchann x 1] values
% nbins: number of interpolation bins per area
% probeType: 'single' (chessboard) or 'dual' (2-shank)

    % --- Step 1: Define shank rules ---
    switch probeType
        case 'NPX2'
            % known x ranges for each shank
            shankRanges = [8 40; 258 290; 508 540; 758 790];
        case 'NPX1'
            % chessboard: only 1 shank, take all channels
            xmin = min(ARRAY1(:,1)); xmax = max(ARRAY1(:,1));
            shankRanges = [xmin xmax];
        otherwise
            error('probeType must be ''single'' or ''dual''');
    end

    nshanks = size(shankRanges,1);

    wrappedVals = cell(nshanks,1);
    wrappedAreas = cell(nshanks,1);
    wrappedY = cell(nshanks,1);
    borders = cell(nshanks,1);

    for s = 1:nshanks
        xlow = shankRanges(s,1);
        xhigh = shankRanges(s,2);

        % --- Step 2: Extract channels for this shank ---
        idxShank = ARRAY1(:,1) >= xlow & ARRAY1(:,1) <= xhigh;
        coords = ARRAY1(idxShank,:);
        areas = CELLARRAY2(idxShank);
        vals = ARRAY3(idxShank);

        if isempty(coords), continue; end

        % --- Step 3: Sort by y position ---
        [~,sortIdx] = sort(coords(:,2),'ascend');
        coords = coords(sortIdx,:);
        areas = areas(sortIdx);
        vals = vals(sortIdx);

        % --- Step 4: Detect area borders ---
        [uniqueAreas,~,ia] = unique(areas,'stable');
        nAreas = numel(uniqueAreas);

        wVals = [];
        wAreas = {};
        wY = [];
        bordersShank = zeros(nAreas,1);

        % --- Step 5: Interpolate within each area ---
        for a = 1:nAreas
            idxArea = ia == a;
            vArea = vals(idxArea);

            % rescale into fixed interval [a, a+1]
            yTarget = linspace(a, a+1, nbins);
            vInterp = interp1(linspace(0,1,numel(vArea)), ...
                               vArea, linspace(0,1,nbins), 'linear');

            wVals = [wVals, vInterp];
            wAreas = [wAreas, repmat(uniqueAreas(a),1,nbins)];
            wY = [wY, yTarget];

            bordersShank(a) = a; % mark border start
        end
        bordersShank(end+1) = nAreas+1; % final border

        wrappedVals{s} = wVals;
        wrappedAreas{s} = wAreas;
        wrappedY{s} = wY;
        borders{s} = bordersShank;
    end
end
