function [str, fun, span, fittedvar] = I_Rs_n_schot(x,y,box,T,varargin)
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
span = [0 max(y)];
V = x(dvec);
I = y(dvec);
I2 = I(1:end-1) + diff(I)/2;
[fitt, gof] = fit(I2,diff(V)./diff(log(I)),'poly1');
n = fitt.p2/(k*T);
Rs = fitt.p1;
fittedvar = [Rs n];
str = ['Fit Range: x =['  num2str(box(1:2)) ']  y=[' num2str(box(3:4)) ']'  newline];
str = [str 'Rs=' num2str(fitt.p1,2)];
str = [str ' \eta=' num2str(n,2)];
str = [str ' R^{2}=' num2str(gof.rsquare,2)];
if ~isempty(varargin)
    str = fittedvar(varargin{:});
end
fun = @(I) Rs*I + n*k*T;
end

