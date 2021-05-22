function [str, fun, span, poly] = I_I0_n_schot_2(x,y,box,R,T,varargin)
if isempty(T)
    T = 296;
end
k = 8.617e-5;
if length(box)==2
    box = [box -inf inf];
elseif isempty(box)
    box = [-inf inf -inf inf];
end
dvec = ~excludedata(x,y,'box',box);
yR = y*R;
V = x(dvec);
I = y(dvec);
[fitt, gof] = fit(V-I*R,log(I./(1-exp(-V/(k*T)))),'poly1');
[~, a] = min(abs(x-box(1)));
[~, b] = min(abs(x-box(2)));
span = [min(0,box(1)-yR(a)) max(max(x-yR),box(2)-yR(b))];
str = ['Fit Range: x =['  num2str([box(1)-yR(a) box(2)-yR(b)]) ']  y=[' num2str(box(3:4)) ']'  newline];
str = [str '\eta=' num2str(1/(k*T*fitt.p1),2)];
str = [str ' I0=' num2str(exp(fitt.p2),2)];
str = [str ' R^{2}=' num2str(gof.rsquare,2)];
n = 1/(k*T*fitt.p1);
lnI0 = fitt.p2;
poly = [exp(fitt.p2) n];
if ~isempty(varargin)
    str = poly(varargin{:});
end
fun = @(V) 1/(k*T*n)*V + lnI0;
end

