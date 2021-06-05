function [StrCells, HandleCells]= PackageNames(varargin)
StrCells = {};
HandleCells = {};
for c = varargin
    ClassNames = {meta.package.fromName(c{:}).ClassList.Name};
    for cc = ClassNames
        in = cc{:};
        StrCells = [StrCells in(find(in=='.',1,'last')+1:end)];
        HandleCells = [HandleCells {eval(['@' in])}];
    end
[StrCells, Indices] = unique(StrCells);
HandleCells = HandleCells(Indices);
end