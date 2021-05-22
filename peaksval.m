function [peakval] = peaksval(s, granularity, cri, m)
if granurality >0
    
    s((end +1 -mod(length(s),granularity)):end,:)=[];
    s  = horzcat(s, (1:length(s))');
    i_smin = min(find(s(:,2)==min(s(:,2))));
    peakval = [];
    i=1;
    sg = mean(vec2mat(s(:,2), granularity), 2);
    
    while 1
    sv = floor(granularity/2):granularity:length(sg);
    sg = horzcat(s(sv,1),sg, sv);
    i_smin = min(find(sg(:,2)==min(sg(:,2))));
    peaksval = sort([peaksval(sg(i_sgmin:-1:1,:), 0, cri, m-1) peaksval(sg((i_sgmin+1):end,:), 0, cri, m-1)

else

end
    
while 1
    if abs(sg(,2)-sg(:,2)

smin = min(sg(:,2));
if s(i,2)-smin>=2*cri
    peakdiv(end+1) = i;
    for k=(i+1):length(s)
        if sg(i,2)-smin<2*cri
                i = k;
                peakdiv(end+1) = k;
                break
        end
    end
end
i=i+1;

while i<length(sg)
    if sg(i,2)-sg(i-1,2)>=cri
        peakdiv(end+1) = i;
        for k=(i+1):length(sg)
            if sg(k,2)-sg(i-1,2)<cri
                i = k;
                peakdiv(end+1) = k;
                break
            end
        end
    end
    i=i;
end

