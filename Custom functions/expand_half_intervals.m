function expanded_regressor = expand_half_intervals(regressor_matrix)
% EXPAND_HALF_INTERVALS
% For each row, find the *central* bout (around middle bin), then extend it
% halfway toward the nearest other bouts (within the same peri-bout window).
% All other values are forced to 0.

[n_bouts, n_bins] = size(regressor_matrix);
expanded_regressor = zeros(size(regressor_matrix));
center_bin = ceil(n_bins/2)+1;

for i = 1:n_bouts
    row = regressor_matrix(i,:)==1;
   

    % --- Find contiguous blocks of 1s ---
    diff_idx = diff([0 row 0]);
    block_starts = find(diff_idx == 1);
    block_ends   = find(diff_idx == -1) - 1;
    if isempty(block_starts)
        block_starts = bout_idx(1);
        block_ends = bout_idx(end);
    end

    % --- Identify the block that includes the center bin ---
    central_block_idx = find(block_starts <= center_bin & block_ends >= center_bin, 1);
    if isempty(central_block_idx)
        % No block overlaps center; pick the one closest to center
        [~,central_block_idx] = min(abs(block_starts - center_bin));
    end
    bout_start = block_starts(central_block_idx);
    bout_end   = block_ends(central_block_idx);

    % --- Find previous and next bouts (within the row) ---
    prev_candidates = block_ends(block_ends < bout_start);
    if isempty(prev_candidates)
        prev_bout_end = 1; % fallback to start of array
    else
        prev_bout_end = max(prev_candidates);
    end
    next_candidates = block_starts(block_starts > bout_end);
    if isempty(next_candidates)
        next_bout_start = n_bins; % fallback to end of array
    else
        next_bout_start = min(next_candidates);
    end

 
    if prev_bout_end==1
        left_extend = prev_bout_end;
    else
        left_extend = round((bout_start + prev_bout_end)/2);
    end

    % --- Compute halfway points for extension ---
  
    if next_bout_start==n_bins
        right_extend = next_bout_start;
    else
        right_extend = round((bout_end + next_bout_start)/2);
    end
    % Right extension: only if a real next bout exists
    
    % --- Fill new row ---
    new_row = zeros(1, n_bins);
    new_row(left_extend:right_extend) = 1;
    expanded_regressor(i,:) = new_row;
end
end