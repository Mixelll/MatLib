function [mat, ind_mat] = shape_in_mat(mat, object, pos, cent)



object_size = size(object);
object_ndims = length(object_size);
object_center = object_size/2 + 0.5*ones(1, object_ndims);


mat_size = size(mat);
mat_size = mat_size(1:object_ndims);
mat_center = mat_size/2 + 0.5*ones(1, ndims(mat));

parity = (rem(mat_size,2)==rem(object_size,2));
if exist('cent', 'var')
     cent = [cent repmat(cent(end),1, object_ndims - length(cent))];
     object_center = object_center .* cent;
else
    cent = ones(1, object_ndims);
end
parity = parity.*~((rem(mat_size,2)==1).*(rem(object_size,2)==1).*(~cent));
parity = parity + (object_size==1) .* (rem(mat_size,2)==0).* (~cent);
mat_center = mat_center + parity;

ind_mat = false(mat_size);
mat_center = num2cell(mat_center);
object_center = num2cell(object_center);
object_size = cellfun(@(x) 0:x-1, num2cell(object_size), 'UniformOutput',false);
pos = num2cell(pos);

% if exist('cent', 'var')
    indices = cellfun(@(a,b,c,d) ceil(a+b-d)*(d>0) + floor(a+b-d)*(d==0) +c, mat_center, pos, object_size, object_center, 'UniformOutput',false);
% else
%     indices = cellfun(@(a,b,c,d) ceil(a+b-d) +c, mat_center, pos, object_size, object_center, 'UniformOutput',false);
% end

if any((mat_size - cellfun(@(x) max(x), indices))<0)
    error('Position of object in matrix exceeds matrix dimensions')
end

indices = [indices repmat({':'}, 1, ndims(mat)-object_ndims)];



ind_mat(indices{:}) = true;
mat(indices{:}) = object;

     

end

