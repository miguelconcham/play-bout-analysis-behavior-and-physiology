function [structure_file, names_file] = find_latest_merged_coherence_files(saving_folder)
%FIND_LATEST_MERGED_COHERENCE_FILES  Newest merged / all-together coherence mat.
%
% After the dual-animal mismatch stitch, the combined pli_struct is saved
% as coherence_structure_merged.mat (or a similar *merged* / *all_together*
% name). This returns that newest structure file and its animal_names
% partner. Satellite pieces (mpfc, B1D11007, NPX2) and the original
% unmerged coherence_structure.mat are ignored.

if nargin < 1 || isempty(saving_folder)
    error('find_latest_merged_coherence_files:Folder', 'saving_folder is required.');
end

listing = [ ...
    dir(fullfile(saving_folder, '*coherence*merged*.mat')); ...
    dir(fullfile(saving_folder, '*coherence*all*together*.mat')); ...
    dir(fullfile(saving_folder, '*coherence*combined*.mat')); ...
    dir(fullfile(saving_folder, '*pli_struct*merged*.mat'))];

keep = false(numel(listing), 1);
for k = 1:numel(listing)
    keep(k) = local_is_merged_structure(listing(k).name);
end
listing = listing(keep);
if isempty(listing)
    error('find_latest_merged_coherence_files:MissingFile', ...
        ['No merged / all-together coherence file in %s. ', ...
        'Run Estimate_coherence_maps.m through the stitch and save ', ...
        'coherence_structure_merged.mat.'], saving_folder);
end

[~, order] = sort([listing.datenum], 'descend');
structure_file = fullfile(listing(order(1)).folder, listing(order(1)).name);
disp(['Using newest merged coherence file: ', structure_file])

names_file = local_partner_names_file(structure_file);
if exist(names_file, 'file') ~= 2
    names_listing = [ ...
        dir(fullfile(saving_folder, '*animal_names*merged*.mat')); ...
        dir(fullfile(saving_folder, '*animal_names*all*together*.mat')); ...
        dir(fullfile(saving_folder, '*animal_names*combined*.mat'))];
    if ~isempty(names_listing)
        [~, names_order] = sort([names_listing.datenum], 'descend');
        names_file = fullfile(names_listing(names_order(1)).folder, ...
            names_listing(names_order(1)).name);
    else
        names_file = '';
    end
end
end

function tf = local_is_merged_structure(filename)
name = lower(filename);
if ~endsWith(name, '.mat')
    tf = false;
    return
end
if contains(name, 'animal_names')
    tf = false;
    return
end
if contains(name, 'mpfc') || contains(name, 'b1d11007') || contains(name, 'npx2')
    tf = false;
    return
end
if ~(contains(name, 'coherence') || contains(name, 'pli_struct'))
    tf = false;
    return
end
tf = contains(name, 'merged') || contains(name, 'all_together') || ...
    contains(name, 'alltogeth') || contains(name, 'all together') || ...
    contains(name, 'combined');
end

function names_file = local_partner_names_file(structure_file)
[folder, stem, ext] = fileparts(structure_file);
if contains(stem, 'coherence_structure')
    names_stem = strrep(stem, 'coherence_structure', ...
        'coherence_structure_animal_names');
else
    names_stem = [stem, '_animal_names'];
end
names_file = fullfile(folder, [names_stem, ext]);
end
