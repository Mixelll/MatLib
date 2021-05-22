function varargout = mat_dim(varargin)

varargout = {};
if nargin<2
    sprintf("returns output of functions accepting repmatted matrices. Then returns the matrices. Matrices are repmatted across sets by index in set"...
        + "\nmat_dim(matSet1,'', matSet2,'', matSetN, range and func pairs, logical([isVec])="...
        + "\nmat_dim(matA1,..., matAn1,'', matB1,...,matBn2,'',..., matNnk, fun1,...,funj,"...
        + "\n..., {range_1A, ...,range_1N},..., {range_jA, ...,range_jN}, logical([isVec]))"...
        + "\nDefinitions:"...
        + "\nSet: DOUBLE/LOGICAL matX1,..., matXn"...
        + "\nfun: FUNCTION_HANDLE function that accepts N-dim matrices"...
        + "\nrange: CELL OF  INT VECTORS. range_i_X of func_i on set X : [indices of matrices from set], e.g. [1:n] or [1 5 7]. Then range_i={range_i_X, range_i_Y}"...
        + "\nIsVec: LOGICAL if set contains 1xN or Nx1 matrices, true if to treat as vec. [set A isvec, set B ismat, sets >=C isvec]==[1 0 1]"...
        + "\nArgument order input does not matter*, they are categorized individually by type."...
        + "\nObey: mat -> double, set separator -> char, fun -> function handler, range -> range vectors in a cell, isVec -> logical"...
        + "\nNo quantity has to be equal: can input N sets, different # of mats in set. Any number of functions. Any number of range cells"...
        + "\nMatrices are expanded across sets, by location in set. set1(mat_k) exapanded by set2(mat_k),..., setN(mat_k)" ...
        + "\nIf set X has n_mat_x>n_mat_y in set Y. Matrices in X with n>ny are expanded by mat_ny in Y"...
        + "\nIf length(isVec)<N sets: isvec for set n>length(isVec) determined by isVec(end). Noninput of isVec -> isVec=true"...
        + "\n*Order between functions and range cells does matter. No cell for fun_>n? fun_>n take all mats from all sets as input")
else
cellcell = {{}};

func = {};
range = {};
set = 1;
for in=varargin
    inn = in{:};
    if isa(inn, 'function_handle')
        func{end+1} = inn;
    elseif isa(inn, 'logical')&&(size(inn,1)==1 || size(inn,2)==1)
        isVec = inn;
    elseif isa(inn, 'double')||isa(inn, 'logical')&&size(inn,1)~=1&&size(inn,2)~=1
        if length(cellcell)<set
            cellcell{end+1} = in;
        else
            cellcell{set}(end+1) = in;
        end
    elseif isa(inn, 'cell')
        range{end+1} = inn;
    elseif isa(inn, 'char')
        set = set +1;
    end
end

% ndimsA = ndims(A_cell{i});
% sizeA = size(A_cell{i});
LF = length(func);
LR = length(range);
LX = cellfun( @length, cellcell); % number of elements in each set
LLX = length(LX); % number of sets
maxLX = max(LX); % # elements in largest set

cell2D = repmat({[]}, maxLX, LLX); % 
cell2D_shifted = cell2D;
range_ind_mat = cell(1, LR);
% ccsize = cell(size(cellcell));
% ccndims = cell(size(cellcell));
c2Dsize = repmat({[]}, maxLX, LLX);
c2Dndims = repmat({[]}, maxLX, LLX);

if exist('isVec', 'var')
    isVec = [isVec repmat(isVec(end),1, LLX - length(isVec))];
else
    isVec = true(1, LLX);
end





 
    for j = 1:length(cellcell)
        jcell = cellcell{j};
        if isVec(j)
            
            for i=1:LX(j)
                if ismatrix(jcell{i})&&(numel(jcell{i})==1)
                    if i==LX(j)&&LX(j)<maxLX
                        cell2D(i:end,j) = jcell(i);
                        c2Dsize(i:end,j) = {[]};
                        c2Dndims(i:end,j) = {[]};
                    else
                        cell2D(i,j) = jcell(i);
                        c2Dsize(i,j) = {[]};
                        c2Dndims(i,j) = {[]};
                    end
                elseif ismatrix(jcell{i})&&(size(jcell{i},1)==1)
                    if i==LX(j)&&LX(j)<maxLX
                        cell2D(i:end,j) = {jcell{i}.'};
                        c2Dsize(i:end,j) = {length(jcell{i})};
                        c2Dndims(i:end,j) = {1};
                    else
                        cell2D(i,j) = {jcell{i}.'};
                        c2Dsize(i,j) = {length(jcell{i})};
                        c2Dndims(i,j) = {1};
                    end
                elseif ismatrix(jcell{i})&&(size(jcell{i},2)==1)
                    if i==LX(j)&&LX(j)<maxLX
                        cell2D(i:end,j) = jcell(i);
                        c2Dsize(i:end,j) = {length(jcell{i})};
                        c2Dndims(i:end,j) = {1};
                    else
                        cell2D(i,j) = jcell(i);
                        c2Dsize(i,j) = {length(jcell{i})};
                        c2Dndims(i,j) = {1};
                    end
                else
                    if i==LX(j)&&LX(j)<maxLX
                        cell2D(i:end,j) = jcell(i);
                        c2Dsize(i:end,j) = {size(jcell{i})};
                        c2Dndims(i:end,j) = {ndims(jcell{i})};
                    else
                        cell2D(i,j) = jcell(i);
                        c2Dsize(i,j) = {size(jcell{i})};
                        c2Dndims(i,j) = {ndims(jcell{i})};
                    end 
                end
            end
                
        else
            for i=1:LX(j)
                if i==LX(j)&&LX(j)<maxLX
                    cell2D(i:end,j) = jcell(i);
                    c2Dsize(i:end,j) = {size(jcell{i})};
                    c2Dndims(i:end,j) = {ndims(jcell{i})};
                else
                    cell2D(i,j) = jcell(i);
                    c2Dsize(i,j) = {size(jcell{i})};
                    c2Dndims(i,j) = {ndims(jcell{i})};
                end
            end
        end
    end

    def_ind_mat = false(maxLX, LLX);
    for j =1:LLX
    def_ind_mat(1:LX(j),j) = true;
    end
    for f=1:LR
        temp_ind_mat = def_ind_mat;
        for j=1:length(range{f})
            temp_ind_mat(:,j) = false;
            temp_ind_mat(range{f}{j},j) = true;
        end
        range_ind_mat{f} = temp_ind_mat;
    end
            
    

for j = 1:length(cellcell)
    for i=1:maxLX
        repvec = [c2Dsize{i,1:j-1} ones(1, c2Dndims{i,j}) c2Dsize{i,j+1:end}];
        if length(repvec)==1
            repvec = [repvec 1];
        end
        cell2D_shifted{i, j} = repmat(shiftdim(cell2D{i,j}, -sum([c2Dndims{i,1:j-1}])), repvec);
    end
end


for f=1:LF
    fout_n = abs(nargout(func{f}));
    if f<=LR

        varargout{f:f-1+fout_n} = func{f}(cell2D_shifted{range_ind_mat{f}}); 
    else
        varargout{f:f-1+fout_n} = func{f}(cell2D_shifted{def_ind_mat});
    end
end

shifted = cell2D_shifted(def_ind_mat);
if size(shifted,1)>size(shifted,2)
    varargout(end+1:end+ sum(def_ind_mat, 'all')) = shifted.';
else
    varargout(end+1:end+ sum(def_ind_mat, 'all')) = shifted;
end
    

end
end

