function [structure_file, names_file] = find_latest_merged_coherence_files(saving_folder, kind)
%FIND_LATEST_MERGED_COHERENCE_FILES  Newest combined coherence mat.
%
% kind:
%   'pairwise' (default) — full GENERATE_COHERNECE_MAPS_STRUCT maps.
%       Newest *merged* / *all_together* / *combined* file, or the
%       all-animal 2026 save (coherence_pli_struct_2026.mat). Ignores
%       satellite pieces (mpfc, B1D11007, NPX2), animal_names, and
%       *lpag_ref* files.
%   'lpag_ref' — GENERATE_COHERENCE_MAPS_LPAG_REF vector saves
%       (coherence_pli_struct_lpag_ref_2026.mat).

if nargin < 1 || isempty(saving_folder)
    error('find_latest_merged_coherence_files:Folder', 'saving_folder is required.');
end
if nargin < 2 || isempty(kind)
    kind = 'pairwise';
end
kind = lower(char(kind));
if ~ismember(kind, {'pairwise', 'lpag_ref'})
    error('find_latest_merged_coherence_files:Kind', ...
        'kind must be ''pairwise'' or ''lpag_ref''.');
end

if strcmp(kind, 'lpag_ref')
    listing = dir(fullfile(saving_folder, '*lpag_ref*.mat'));
else
    listing = [ ...
        dir(fullfile(saving_folder, '*coherence*merged*.mat')); ...
        dir(fullfile(saving_folder, '*coherence*all*together*.mat')); ...
        dir(fullfile(saving_folder, '*coherence*combined*.mat')); ...
        dir(fullfile(saving_folder, '*pli_struct*merged*.mat')); ...
        dir(fullfile(saving_folder, '*pli_struct*2026*.mat')); ...
        dir(fullfile(saving_folder, '*coherence*2026*.mat'))];
end

if ~isempty(listing)
    [~, unique_idx] = unique({listing.name}, 'stable');
    listing = listing(unique_idx);
end

keep = false(numel(listing), 1);
for k = 1:numel(listing)
    keep(k) = local_is_structure(listing(k).name, kind);
end
listing = listing(keep);
if isempty(listing)
    if strcmp(kind, 'lpag_ref')
        error('find_latest_merged_coherence_files:MissingFile', ...
            ['No LPAG-ref coherence file in %s. ', ...
            'Run Estimate_coherence_maps_lpag_ref.m and save ', ...
            'coherence_pli_struct_lpag_ref_2026.mat.'], saving_folder);
    else
        error('find_latest_merged_coherence_files:MissingFile', ...
            ['No merged / all-together / 2026 coherence file in %s. ', ...
            'Run Estimate_coherence_maps.m and save ', ...
            'coherence_pli_struct_2026.mat.'], saving_folder);
    end
end

[~, order] = sort([listing.datenum], 'descend');
structure_file = fullfile(listing(order(1)).folder, listing(order(1)).name);
disp(['Using newest ', kind, ' coherence file: ', structure_file])

names_file = local_partner_names_file(structure_file);
if exist(names_file, 'file') ~= 2
    names_listing = dir(fullfile(saving_folder, '*animal_names*.mat'));
    keep_names = false(numel(names_listing), 1);
    for k = 1:numel(names_listing)
        keep_names(k) = local_is_names_fallback(names_listing(k).name, kind);
    end
    names_listing = names_listing(keep_names);
    if ~isempty(names_listing)
        [~, names_order] = sort([names_listing.datenum], 'descend');
        names_file = fullfile(names_listing(names_order(1)).folder, ...
            names_listing(names_order(1)).name);
    else
        names_file = '';
    end
end
end

function tf = local_is_structure(filename, kind)
name = lower(filename);
if ~endsWith(name, '.mat')
    tf = false;
    return
end
if contains(name, 'animal_names')
    tf = false;
    return
end
if ~(contains(name, 'coherence') || contains(name, 'pli_struct'))
    tf = false;
    return
end
is_lpag_ref = contains(name, 'lpag_ref');
if strcmp(kind, 'lpag_ref')
    tf = is_lpag_ref;
    return
end
if is_lpag_ref
    tf = false;
    return
end
if contains(name, 'mpfc') || contains(name, 'b1d11007') || contains(name, 'npx2')
    tf = false;
    return
end
tf = contains(name, 'merged') || contains(name, 'all_together') || ...
    contains(name, 'alltogeth') || contains(name, 'all together') || ...
    contains(name, 'combined') || ...
    (contains(name, '2026') && contains(name, 'pli_struct'));
end

function tf = local_is_names_fallback(filename, kind)
name = lower(filename);
if ~endsWith(name, '.mat') || ~contains(name, 'animal_names')
    tf = false;
    return
end
is_lpag_ref = contains(name, 'lpag_ref');
if strcmp(kind, 'lpag_ref')
    tf = is_lpag_ref;
    return
end
if is_lpag_ref
    tf = false;
    return
end
tf = contains(name, 'merged') || contains(name, 'all_together') || ...
    contains(name, 'alltogeth') || contains(name, 'all together') || ...
    contains(name, 'combined') || contains(name, '2026');
end

function names_file = local_partner_names_file(structure_file)
[folder, stem, ext] = fileparts(structure_file);
if strcmp(stem, 'coherence_pli_struct_2026')
    names_file = fullfile(folder, ['coherence_pli_structure_animal_names2026', ext]);
    return
end
if strcmp(stem, 'coherence_pli_struct_lpag_ref_2026')
    names_file = fullfile(folder, ['coherence_pli_struct_lpag_ref_animal_names_2026', ext]);
    return
end
if contains(stem, 'coherence_structure')
    names_stem = strrep(stem, 'coherence_structure', ...
        'coherence_structure_animal_names');
else
    names_stem = [stem, '_animal_names'];
end
names_file = fullfile(folder, [names_stem, ext]);
end
