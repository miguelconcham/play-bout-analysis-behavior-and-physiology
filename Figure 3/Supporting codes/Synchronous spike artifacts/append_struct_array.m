function collection = append_struct_array(collection, items)
%APPEND_STRUCT_ARRAY  Append structs that do not share every field.
%
% Older session mats are missing fields added later (pairs, coverage,
% nonsync_mean). orderfields errors when the field lists differ in length.
% Missing fields are filled with [] so the mats can still be combined.

if isempty(items)
    return
end
for index = 1:numel(items)
    collection = append_one(collection, items(index));
end
end

function collection = append_one(collection, item)
if isempty(collection)
    collection = item;
    return
end
names = union(fieldnames(collection), fieldnames(item), 'stable');
collection = ensure_fields(collection, names);
item = ensure_fields(item, names);
collection(end + 1) = orderfields(item, collection); %#ok<AGROW>
end

function value = ensure_fields(value, names)
for index = 1:numel(names)
    if ~isfield(value, names{index})
        [value.(names{index})] = deal([]);
    end
end
value = orderfields(value, names);
end
