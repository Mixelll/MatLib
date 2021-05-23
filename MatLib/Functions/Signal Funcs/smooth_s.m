function sy = smooth_s(sx, sy, r, method, nmins)

allsmindex = [];
if size(sx, 2)==1
    sx = repmat(sx, 1, size(sy, 2));
end
for k=1:size(sy, 2)
    syk = sy(:,k);
    syk = syk(1:sum(~isnan(syk)));
    syksorted = sort(sy(:,k));
    ykmin = mean(syksorted(1:nmins));
    bar = r + ykmin;
    smindex = [];
    in = 0;
    for i=1:length(syk)
        if (syk(i)<bar)&&~in
            smindex(end+1,1) = i;
            in = 1;
        end
        if (syk(i)>=bar)&&in
            smindex(end,2) = i-1;
            in = 0;
        end
    end
    
    if (length(smindex)<2)||(smindex(end,2)==0)
        smindex(end,2) = length(syk);
    end
        for i=1:size(smindex,1)
            if smindex(i,2)>smindex(i,1)
                sy(smindex(i,1):smindex(i,2), k) = method(sx(smindex(i,1):smindex(i,2), k), syk(smindex(i,1):smindex(i,2)));
            end
        end
        
    allsmindex(1:(2*size(smindex, 1)), k) = smindex(:);
end
end

