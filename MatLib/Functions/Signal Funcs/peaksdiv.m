function [peakdiv] = peaksdiv(s, cri)
dsdx = dif(s, 1);
x = dsdx(:,1);
dydx = dsdx(:,2);
peakdiv = x(abs(dydx)>cri);

end

