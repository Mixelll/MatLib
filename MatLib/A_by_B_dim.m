function varargout = A_by_B_dim(varargin)
varargout = {};
A_cell = {};
B_cell = {};
setA = 1;
func_indicator = 0;
for in=varargin
    if isa(in{:}, 'function_handle')
        func = in{:};
        func_indicator = 1;
    elseif isa(in{:}, 'double')&&(length(in{:})==1)&&((in{:}==1)||(in{:}==0))
        isVec = in{:};
    elseif isa(in{:}, 'double')&&setA
        A_cell(end+1) = in;
    elseif isa(in{:}, 'double')
        B_cell(end+1) = in;
    elseif isa(in{:}, 'cell')
        range = in{:}{:};
    else
        setA = 0;
    end
end
  
% ndimsA = ndims(A_cell{i});
% sizeA = size(A_cell{i});
LA = length(A_cell); 
LB = length(B_cell);

if exist('isVec', 'var')&&isVec

    for i=1:LA
        if ismatrix(A_cell{i})&&(size(A_cell{i},2)==1)
            A_cell{i} = A_cell{i}';
        end
    end
    
    for i=1:LB
        if ismatrix(B_cell{i})&&(size(B_cell{i},1)==1)
            B_cell{i} = B_cell{i}';
        end
    end
end

inner_dims = A_cell;
outer_dims = B_cell;
L_inner = length(inner_dims);
L_outer = length(outer_dims);

for i=1:L_inner

    A_cell{i} = repmat(inner_dims{i}, [ones(1, ndims(inner_dims{i})) size(outer_dims{min(i,L_outer)}) ]);
end

for i=1:L_outer
    
    B_cell{i} = repmat(shiftdim(outer_dims{i}, -ndims(inner_dims{min(i,L_inner)})), size(inner_dims{min(i,L_inner)}));
end



if func_indicator
    if exist('range', 'var')
        varargout{1} = func(A_cell{1:range(1)}, B_cell{1:range(end)});
    else
        varargout{1} = func(A_cell{1:end}, B_cell{1:end});
    end
end

varargout = [varargout cell(1, LA+LB)];

if  exist('isVec', 'var')&&isVec
    A_cell = cellfun(@shiftdim, A_cell, 'UniformOutput', false);
    B_cell = cellfun(@shiftdim, B_cell, 'UniformOutput', false);
    if func_indicator
        varargout{1} = shiftdim(varargout{1});
    end
end

varargout(func_indicator+1     : LA    +func_indicator) = A_cell;
varargout(func_indicator+1 +LA : LA+LB +func_indicator) = B_cell;
varargout{end+1} = size(B_cell{1});

end

