function base = analysis_figure_basename(figure_name, figure_number)
%ANALYSIS_FIGURE_BASENAME  File name for one saved analysis figure.
%
% Empty names use the figure number. Characters that are illegal in a
% Windows file name are removed so the same folder can be written on the
% lab machine.

if nargin < 2 || isempty(figure_number)
    figure_number = 1;
end
if nargin < 1 || isempty(figure_name)
    figure_name = sprintf('figure %d', figure_number);
end
base = char(figure_name);
base = strrep(base, ':', ' - ');
base = regexprep(base, '[\\/:*?"<>|]', '-');
base = regexprep(base, '\s+', ' ');
base = strtrim(base);
if isempty(base)
    base = sprintf('figure %d', figure_number);
end
end
