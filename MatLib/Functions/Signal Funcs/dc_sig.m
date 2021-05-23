function sy = dc_sig(sx, sy, r, DCns, varargin)
if ~isempty(varargin)&&isa(varargin{1}, 'numeric')
    user_def_value = varargin{1};
end

if ~exist('user_def_value', 'var')
    nsegments = DCns(1);
    nmins = DCns(2);
    if size(sx, 2)==1
        sx = repmat(sx, 1, size(sy, 2));
    end
    for k=1:size(sy, 2)
        syk = sy(:,k);
        syk = syk(1:sum(~isnan(syk)));
        syksorted = sort(syk);
        ykmin = mean(syksorted(1:nmins));
        bar = r + ykmin;
        lsyk = length(syk);
        lsegment = floor(lsyk/nsegments);
        segindex = horzcat((1:lsegment:((nsegments-1)*lsegment +1))', [lsegment:lsegment:(nsegments-1)*lsegment lsyk]');
        dcindex = [];
        dcvalue = [];
        for i=1:size(segindex,1)
            segmentunsorted = syk(segindex(i,1):segindex(i,2));
            segmentsorted = sort(segmentunsorted);
            segmentsorted(nmins);
            dcseglogical = segmentunsorted<=min(segmentsorted(nmins), bar);
            dcseg = mean(segmentunsorted(segmentunsorted<=min(segmentsorted(nmins), bar)));
            dcsegiarray = find(dcseglogical);
            dcsegi = floor(mean(dcsegiarray));
            if ~isnan(dcsegi)
                dcindex(end+1) = dcsegi + (i-1)*lsegment;
                dcvalue(end+1) = dcseg;
            end

        end
        if length(dcindex)>1
            dcy = interp1(sx(dcindex(:),k), dcvalue, sx(1:lsyk,k), 'linear', 'extrap');
        else
            dcy = ykmin*ones(lsyk, 1);
        end
        sy(1:lsyk,k) = dcy;
    end
else
    sy(~isnan(sy)) = user_def_value;
end
end
