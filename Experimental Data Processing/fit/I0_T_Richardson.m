function [str, fun, span, fittedvar] = I0_T_Richardson(x,y,box,A,varargin)
k = 8.617e-5;
if length(box)==2
    box = [box -inf inf];
elseif isempty(box)
    box = [-inf inf -inf inf];
end

dvec = ~excludedata(x,y,'box',box);
span = [1/max(x) 1/min(x)];
T = x(dvec);
I0 = y(dvec);
[fitt, gof] = fit(1./T,log(I0./T.^2),'poly1');
VB = -fitt.p1*k;
As = exp(fitt.p2)/A;
fittedvar = [As VB];
str = ['Fit Range: x =['  num2str(box(1:2)) ']  y=[' num2str(box(3:4)) ']'  newline];
str = [str 'As=' num2str(As,2)];
str = [str ' VB=' num2str(VB,2)];
str = [str ' R^{2}=' num2str(gof.rsquare,2)];
if ~isempty(varargin)
    str = fittedvar(varargin{:});
end
fun = @(oneoverT) -VB/k.*oneoverT + log(A.*As);
end

