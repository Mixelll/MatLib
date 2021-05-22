function [scells, refi] = data_folder_read(path, sortby, ind, tparam, saveop, varargin)

baseline_search_str = ["ref", "baseline", "DC"];

CharOrString = @(s) ischar(s) || isstring(s);
p = inputParser;
p.KeepUnmatched=true;
p.addParameter('FileType', '', CharOrString);
p.parse(varargin{:});

%max=0;

if isempty(ind)
    index = '';
    trail = '';
elseif ind(1)=='m'
    index = '';
    trail = ind(2:end);
%     max = 1;
else
    index = num2str(ind(1));
    trail = ind(2:end);
end
lind = length(index);
ltrail = length(trail);

if isempty(p.Results.FileType)
    dir_string = strcat(path, '\*', index, trail);
else
    dir_string = strcat(path, '\*', index, trail, '.', p.Results.FileType);
end
dircontents = dir(dir_string); % this is the structure.fiels that contains the ind&type filtered directory files
dircontents = struct2table(dircontents);
dircontents(dircontents.isdir,:) = [];
dircontents = sortrows(dircontents, 'name');
dircontents = table2struct(dircontents);
n = length(dircontents);
maxi = zeros(1,n);

% some sorting. sortby - location of index from end of name 
if isempty(sortby)
    for i=1:n
        cnamep = dircontents(i).name;
        cnamef = dircontents(min(i+1,n)).name;
        dircontents(i).index = str2double(cnamep(end -4 -ltrail));
        dircontents(i).index;
        if ~strcmp(cnamep(1:end -lind -ltrail),cnamef(1:end -lind -ltrail))
            maxi(i) = 1;
        end
    end
    maxi(end) = 1;
    dircontents = dircontents(logical(maxi));
else
    for i=1:n
        cnamep = dircontents(i).name;
        sindex = '';
        for s=sortby
            st = cnamep(end-3-s);
            if (st>='0' && st<='9') || st=='-' || st=='+' 
                sindex = [sindex st];
            end

        end
        dircontents(i).index = str2num(sindex);
    end
    [~,ss] = sort([dircontents.index]);
    dircontents = dircontents(ss);
end




% populate cells with matrixed file contents
n = length(dircontents);
scells = cell(3,n);
refi = 0;
tparamflag = 0;
for i=1:n
    pathi = strcat(path, '\', dircontents(i).name);
    fullname = dircontents(i).name;
    if contains(fullname, baseline_search_str,'IgnoreCase',true)
        refi = i;
    end  
    scells{1,i} = fullname(1:end-4);
    scells{2,i} = readmatrix(pathi, varargin{:});
    for t = tparam
        tparamflag = 1;
        tcells = {};
        st = scells{1,i};
        found = strfind(st, t{:});
        st = st(found(end):end);
        if ~isempty(found)
            numt = (st>='0' & st<='9') | st=='-' | st=='+'  | st=='.';
            tcells = [tcells; t {str2double(st(find(numt,1):find(diff(numt)==-1,1)))}];
        end
    end
    if tparamflag
        scells{3, i} = tcells;
    end
end
% move baseline to end
if refi
    temp = scells(:,refi);
    scells(:,refi) = [];
    scells(:,end+1) = temp;
end

if saveop==1
    save([ path '\saved.mat'], 'scells');
elseif saveop==2
    for i=1:size(scells,2)
        temp = scells{1,i};
        save([ path '\' scells{1,i} '.mat'], 'temp');
    end
end
        
end

