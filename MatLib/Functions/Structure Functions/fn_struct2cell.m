function result_cell = fn_struct2cell(Xname, varargin)
% function fn_structdisp Xname
% function fn_structdisp(X)
%---
% Recursively display the content of a structure and its sub-structures
%
% Input:
% - Xname/X     one can give as argument either the structure to display or
%               or a string (the name in the current workspace of the
%               structure to display)
%
% A few parameters can be adjusted inside the m file to determine when
% arrays and cell should be displayed completely or not

% Thomas Deneux
% Copyright 2005-2012
% Michael Leitsin
% Copyright 2020-2021
p = inputParser;
p.addParameter('ignore', {})
p.addParameter('write', false) % For object handles: Return only writable properties
p.parse(varargin{:})

result_cell = {};
if ~isempty(Xname)
    if ischar(Xname)
        X = evalin('caller',Xname);
    else
        X = Xname;
        Xname = inputname(1);
    end

    if ~(isstruct(X) || isobject(X)), warning('argument should be a structure or the name of a structure or an object'), result_cell={''; X; ''; ''}; return, end
    result_cell = rec_struct2cell(Xname,X,result_cell, p.Results.ignore, p.Results.write);
    for i = 1:size(result_cell,2)
        locs = strfind(result_cell{1,i}, '.');
        result_cell{3,i} = result_cell{1,i}(locs(1):end);
        result_cell{4,i} = result_cell{1,i}(locs(end)+1:end);
    end
end

%---------------------------------
function result_cell = rec_struct2cell(Xname,X,result_cell, ignore, write)
%---

%-- PARAMETERS (Edit this) --%

ARRAYMAXROWS = 100000;
ARRAYMAXCOLS = 100000;
ARRAYMAXELEMS = 1000000;
CELLMAXROWS = 100000;
CELLMAXCOLS = 100000;
CELLMAXELEMS = 1000000;
CELLRECURSIVE = false;

%----- PARAMETERS END -------%


if isstruct(X) || isobject(X)
    if isobject(X) && write
        try 
            F = fieldnames(set(X));
        catch
            try 
                F = fieldnames(X);
            catch
                F = {};
            end
        end
    else
        F = fieldnames(X);
    end
    F = F(cellfun(@(c) ~any(strcmpi(c,ignore)),F));
    nsub = length(F);
    Y = cell(1,nsub);
    subnames = cell(1,nsub);
    for i=1:nsub
        f = F{i};
        Y{i} = X.(f);
        subnames{i} = [Xname '.' f];
    end
elseif CELLRECURSIVE && iscell(X)
    nsub = numel(X);
    s = size(X);
    Y = X(:);
    subnames = cell(1,nsub);
    for i=1:nsub
        inds = s;
        globind = i-1;
        for k=1:length(s)
            inds(k) = 1+mod(globind,s(k));
            globind = floor(globind/s(k));
        end
        subnames{i} = [Xname '{' num2str(inds,'%i,')];
        subnames{i}(end) = '}';
    end
elseif iscell(X)
    s = size(result_cell,2);
    result_cell{1,s+1} = Xname;
    result_cell{2,s+1} = X;
    return
else
    return
end

 for i=1:nsub
    a = Y{i};
    if isstruct(a) || isobject(a) && ~islogical(a) && ~isa(a, 'sym')
        if length(a)==1
            result_cell = rec_struct2cell(subnames{i},a,result_cell, ignore, write);
        else
            for k=1:length(a)
                result_cell = rec_struct2cell([subnames{i} '(' num2str(k) ')'],a(k),result_cell, ignore, write);
            end
        end
    elseif iscell(a)
        if size(a,1)<=CELLMAXROWS && size(a,2)<=CELLMAXCOLS && numel(a)<=CELLMAXELEMS
            result_cell = rec_struct2cell(subnames{i},a,result_cell, ignore, write);
        end
    elseif size(a,1)<=ARRAYMAXROWS && size(a,2)<=ARRAYMAXCOLS && numel(a)<=ARRAYMAXELEMS
        s = size(result_cell,2);
        result_cell{1,s+1} = subnames{i};
        result_cell{2,s+1} = a;
%         disp([subnames{i} ':'])
%         disp(a)
    end
end
