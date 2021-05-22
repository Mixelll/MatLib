function result_cell = s_noise(labels, sx, sy, smooth_init, range, scale, units)
if size(sx, 2)==1
    sx = repmat(sx, 1, size(sy, 2));
end

sysmoothed = smooth_init(sx, sy);
serror = ((sy - sysmoothed)).^2;
serror_byscale = serror/scale;
[~, result_cell] = signal_area(labels, sx, serror_byscale, range, scale, '');
% plot(sx, sy, sx, sysmoothed, sx, serror_byscale)

for i=1:length(range)/2
    rangewidth = range(2*i) - range(2*i -1);
    result_cell(2:end,i+1) = num2cell(cell2mat(result_cell(2:end,i+1))./rangewidth);
end   
    % valuevector = sysmoothed(sysmoothed>0);
% valuevector = sort(valuevector(:));
% sysmoothed_no0(sysmoothed<=0) = mean(valuevector(1:(size(sy, 1)*size(sy, 2)/2)));
% serror_byvalue = ((sy - sysmoothed)./sysmoothed_no0).^2;
% signal_area(labels, sx, serror_byvalue, range, units)
% plot(sx, sy, sx, sysmoothed, sx, serror_byvalue)
% result_array = zeros(size(integral,2), length(range)/2 +1);
% result_cell = cell(size(integral,2)+1, length(range)/2 +1);
% result_cell{1,1} = units;
% result_cell(2:end,1) = labels';
% for k=1:size(syk,2)
%     sxk = sx(left:right, k)
%     syk = syk(1:sum(~isnan(syk)));
%     smy = smooths_init(sxk(:,k), syk(:,k));
%     for i=1:length(range)/2
%         [~, left] = min(abs(sxk(:,k) - range(2*i -1)));
%         [~, right] = min(abs(sxk(:,k) - range(2*i)));
%         vecx = sxk(left:right, k);
%         vecy = syk(left:right, k);
%         vecsmy = smy(left:right, k);
%         vecsmy_no0 = vecsmy;
%         vecsmy_no0(vecsmy_no0<=0) = min(vecsmy(vecsmy>0));
%         serror_byvalue = ((vecy - smy)./vecsmy_no0).^2;
%         vec_integ(vecx, serror_byvalue)
%         result_array(k,i) = vec_integ(vecx, vecy);
%         result_cell{k+1,i+1} = result_array(k,i);
%         result_cell{1,i+1} = [num2str(range(2*i -1)) '-' num2str(range(2*i))];
%     end
% end
    

end

