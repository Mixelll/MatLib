function [OUT, OK, Pairs] = StructrureFieldsMenu(p, varargin)
OK = true;
OUT = struct;
Pairs = {};
Cells = fn_struct2cell(p);
Title = 'Input Parameters';
Prompt = Cells(1,:);
VararginRest = {};
for c = varargin
    if ischar(c{:}) || isstring(c{:})
        Title = c{:};
    elseif  iscell(c{:})
        Prompt = c{:};
    else
        VararginRest(end+1) = c;
    end
end
if numel(VararginRest)<2
    ParseInMenu = @parse_num_cell_sym2char;
    ParseOutMenu = @parse_str2num_cell;
else
    ParseInMenu = VararginRest{1};
    ParseOutMenu = VararginRest{2};
end
opts.Resize = 'off';
opts.WindowStyle = 'normal';
opts.Interpreter = 'none';
answer = inputdlg(Prompt, Title, [1 30+numel(Title)],ParseCells(Cells(2,:),ParseInMenu), opts);
if ~isempty(answer)
    NewValues = ParseCells(answer,ParseOutMenu);
    Cells(2,:) = NewValues;
    for c = Cells
        Temp = c{2};
        eval(['OUT' c{3} '=Temp;'])
    end
    StructCell = fn_struct2cell(OUT);
    Pairs = [StructCell(4,:) ; StructCell(2,:)];
else
    OK = false;
end
end

