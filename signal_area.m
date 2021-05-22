function [result_cell_native, result_cell_scaled, result_cell_scaled_fraction, result_array_native] = signal_area(labels, sx, sy, range, scale, units)

if size(sx, 2)==1
    sx = repmat(sx, 1, size(sy, 2));
end
syinteg = vec_integ(sx, sy);
lsyinteg = length(syinteg);
result_array_native = zeros(size(syinteg,2), length(range)/2);
result_cell_native = cell(size(syinteg,2)+1, length(range)/2 +1);
result_cell_native{1,1} = units;
result_cell_native(2:end,1) = labels';
result_cell_scaled = result_cell_native;
result_cell_scaled{1,1} = [units ' scaled'];
for k=1:size(syinteg,2)
    for i=1:length(range)/2
        [~, left] = min(abs(sx(:,k) - range(2*i -1)));
        [~, right] = min(abs(sx(:,k) - range(2*i)));
        left = min(left, lsyinteg);
        right = min(right, lsyinteg);
        right = right*(right<=length(syinteg)) + length(syinteg)*(right>length(syinteg));
        result_array_native(k,i) = syinteg(right, k) - syinteg(left, k);
        result_cell_native{k+1,i+1} = result_array_native(k,i);
        result_cell_native{1,i+1} = [num2str(range(2*i -1)) '-' num2str(range(2*i))];
        
        result_cell_scaled{k+1,i+1} = result_array_native(k,i)/scale;
        result_cell_scaled{1,i+1} = [num2str(range(2*i -1)) '-' num2str(range(2*i))];
    end
end
result_cell_scaled_fraction = result_cell_scaled;
for j=2:size(result_cell_scaled,2)
    last_val = result_cell_scaled{end, j};
    result_cell_scaled_fraction(2:end, j) = cellfun(@(x) x/last_val,result_cell_scaled(2:end, j),'UniformOutput',false);
end

end

