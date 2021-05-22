function [str, fun, span, poly] = I_I0_n_schot_1(x,y,box,T,varargin)
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
span = [0 max(x)];
V = x(dvec);
I = y(dvec);
[fitt, gof] = fit(V,log(I./(1-exp(-V/(k*T)))),'poly1');
str = ['Fit Range: x =['  num2str(box(1:2)) ']  y=[' num2str(box(3:4)) ']'  newline];
str = [str '\eta=' num2str(1/(k*T*fitt.p1),2)];
str = [str ' I0=' num2str(exp(fitt.p2),2)];
str = [str ' R^{2}=' num2str(gof.rsquare,2)];
n = 1/(k*T*fitt.p1);
LnI0 = fitt.p2;
poly = [fitt.p2 fitt.p1];
if ~isempty(varargin)
    str = poly(varargin{:}+1);
end
fun = @(V) 1/(k*T*n)*V + LnI0;
end

