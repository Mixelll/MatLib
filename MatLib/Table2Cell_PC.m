function CellsOut = Table2Cell_PC(DataTableIn, Param, varargin)

    
VarNames = DataTableIn.Properties.VariableNames;
NameVar = VarNames(contains(VarNames,'name', 'IgnoreCase',true)); NameVar = NameVar{1};
DataVar = VarNames(strcmpi(VarNames,'data')); DataVar = DataVar{1};

if ~isempty(varargin)
    StructVar = varargin{1};
    StructVar = VarNames(strcmpi(VarNames,StructVar)); StructVar = StructVar{1};
    StructParams = unique(DataTableIn.(StructVar)).';
    CellsOut = struct;
else
    StructParams = 1;
end
ParamFlag = false;
if ~isempty(Param)
    if isstring(Param)
        ParamFlag = true;
        Param = cellstr(Param);
    elseif ~iscell(Param)
        ParamFlag = true;
        Param = {Param};
    elseif ~isempty(Param{1})
        ParamFlag = true;
    end
    Param = Param(~cellfun(@isempty, Param));
end
    
DataTable = DataTableIn;
for SP = StructParams
    if ~isnumeric(SP)
        DataTable = DataTableIn(DataTableIn.(StructVar)==SP,:);
    end
    Cells = DataTable.(NameVar).';
    if ~iscell(Cells)
        Cells = num2cell(Cells);
    end
    for i = 1:size(Cells,2)
        Cells{i} = convertStringsToChars(Cells{i});
    end
    Data = DataTable.(DataVar).';
    if ~iscell(Data)
        Data = num2cell(Data);
    end
    Cells = [Cells ; Data];
    if ParamFlag
        ParamCellAtot = cell(length(Param), size(Cells,2));
        for i = 1:length(Param) 
            ParamCellA = DataTable.(Param{i}).';
            if ~iscell(ParamCellA)
                ParamCellA = num2cell(ParamCellA);
            end
            ParamCellAtot(i,:) = ParamCellA;
        end
        ParamCellB = cell(1, size(Cells,2)); 
        for i = 1:size(Cells,2)
            ParamCellB{i} = [Param' ParamCellA(:,i)];
        end
        Cells = [Cells ; ParamCellB];
    end
    if ~isnumeric(SP)
        CellsOut.(SP) = Cells;
    else
        CellsOut = Cells;
    end
end

