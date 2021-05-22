function [pos] = peakpos(s, granularity_x, granularity_y)
% s((end +1 -mod(length(s),granularity_x)):end,:)=[];
s  = horzcat(s, (1:length(s))');
wall = zeros(length(s), 1);
smin = min(s(:,2));
smax = max(s(:,2));
d = (smax - smin)/granularity_y;
ds = d;
sp=s;
shelf = ones(granularity_x,1);
fail = 0;
while fail<20&&(ds + smin)<1200
    wall_l = wall;
    for i=1:length(sp)
%         i
%         ((length(sp) - i)>=granularity_x)
%         all(((ds + smin)*shelf-sp(i:i+(granularity_x - 1),2))>0)
%         all(wall(i:(i+granularity_x - 1))==0)
        %if ((length(sp) - i)>=granularity_x)
        if ((length(sp) - i)>=granularity_x)&&all(((ds + smin)*shelf-sp(i:i+(granularity_x - 1),2))>0)&&all(wall(i:(i+granularity_x - 1))==0)
            wall(i:(i+(granularity_x+1)-1)) = ds + smin;
        %end
        end
    end
    if all(wall_l==wall)
        fail = fail + 1;
    else
        fail = 0;
    end
    ds = ds + d;
end
pos  = horzcat(s, wall);
