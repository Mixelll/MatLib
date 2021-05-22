function [peakdiv] = peaks(s, cri)
dsdx = dif(s, 1);
x = dsdx(:,1);
dydx = dsdx(:,2);
peak = x(abs(dydx)>cri);

end

