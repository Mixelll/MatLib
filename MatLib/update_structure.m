function [original_struct] = update_structure(original_struct,updated_struct,varargin)

p = inputParser;
p.addParameter('ignore', {}) % field branches to ignore
p.addParameter('write', true) % For object handles: Return only writable properties
p.addParameter('new', false) % True for add new properies. False to update only existing properties
p.addParameter('onlynew', false) % Only add new properties, don't update existing
p.addParameter('case', false) % compare properties case sensivity
p.parse(varargin{:})
new = p.Results.new;
if p.Results.onlynew
    new = true;
end

original_cell = fn_struct2cell(original_struct, 'ignore',p.Results.ignore, 'write',p.Results.write);
updated_cell = fn_struct2cell(updated_struct, 'ignore',p.Results.ignore, 'write',p.Results.write);
if ~isempty(fieldnames(updated_struct))
    if new
        for c = updated_cell
            if ~isempty(original_cell)
                if p.Results.case
                    field_cmp = strcmp(original_cell(3,:), c{3});
                else
                    field_cmp = strcmpi(original_cell(3,:), c{3});
                end
            end
            if p.Results.onlynew
                if isempty(original_cell) || all(~field_cmp)
                    eval(['original_struct' c{3} '=' c{1} ';'])
                end
            else
                if ~isempty(original_cell) && any(field_cmp)
                    eval(['original_struct' original_cell{3,field_cmp} '=' c{1} ';'])
                else
                    eval(['original_struct' c{3} '=' c{1} ';'])
                end
            end
        end
    else
        for c = original_cell
            if p.Results.case
            	field_cmp = strcmp(updated_cell(3,:), c{3});
            else
                field_cmp = strcmpi(updated_cell(3,:), c{3});
            end
            if any(field_cmp)
                eval(['original_struct' c{3} '=' updated_cell{1,field_cmp} ';'])
            end
        end
    end
end
end 

