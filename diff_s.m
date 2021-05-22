function [s] = diff_s(s, n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
for i=1:n
    d = diff(s);
    dydx = d(:,2)./d(:,1);
    s = horzcat(s(1+(1+(-1)^i)/2:end-(1+(-1)^(i+1))/2,1), dydx);
end
end

