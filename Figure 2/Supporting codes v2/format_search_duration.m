function text = format_search_duration(seconds)
%FORMAT_SEARCH_DURATION  Short clock text for progress lines.

if ~isfinite(seconds)
    text = 'unknown';
    return
end
seconds = max(0, seconds);
if seconds < 90
    text = sprintf('%.0f s', seconds);
elseif seconds < 3600
    text = sprintf('%.1f min', seconds / 60);
else
    text = sprintf('%.1f h', seconds / 3600);
end
end
