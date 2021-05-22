function [out] = stat_analysis(data)
title = {'mean' 'std' 'skewness' 'kurtosis'};
weighted = repmat(data(:,1), 1, size(data,2)-1).*data(:,2:end)./repmat(mean(data(:,2:end)), size(data,1), size(data,2)-1);

std_arr = zeros(1, size(data,2)-1);
for i = 2:size(data,2)
    std_arr(1,i-1) = std(data(:,1), data(:,i));
end
out = {title{:}; mean(weighted) std_arr skewness(weighted) kurtosis(weighted)};

end

