function [str, fun, span, FittedParamCell] = str_poly_fit(x,y,n,varargin)
p = inputParser;
p.addParameter('range', [-inf inf -inf inf], @isnumeric); % fitting range
p.addParameter('neglect', 0, @isnumeric); % neglect coefficients that contribute less than specified fraction
p.addParameter('params', {}, @iscell); % parameter aliases to put in FittedParamCell{3,:}
p.parse(varargin{:});

if length(p.Results.range)==2
    range = [p.Results.range -inf inf];
else
    range = p.Results.range;
end

dvec = ~excludedata(x,y,'box',range) & ~isnan(y);
x = x(dvec);
y = y(dvec);
xMid = round(length(x)/2);

str = [];
[PolyCoef, S] = polyfit(x,y,n);
R2 = 1 - (S.normr/norm(y - mean(y)))^2;
PowerVec = n:-1:0;

p_keep_logical = abs(PolyCoef.*x(xMid).^PowerVec) > p.Results.neglect*max(abs(PolyCoef.*x(xMid).^PowerVec));
PolyCoef = PolyCoef(p_keep_logical);
PowerVec = PowerVec(p_keep_logical);
PolyCoef = PolyCoef(end:-1:1);
PowerVec = PowerVec(end:-1:1);

for i=1:length(PowerVec)
    if PowerVec(i)>1
        str = [str num2str(PolyCoef(i)) '*x^{' num2str(PowerVec(i)) '}'];
    elseif PowerVec(i)==1
        str = [str num2str(PolyCoef(i)) '*x'];
    else
        str = [str num2str(PolyCoef(i))];
    end
    if i<length(PowerVec)
        str = [str ' +'];
    end
end
str = [str newline 'R^{2} = ' num2str(R2,3)];
FittedParamCell = cell(3,length(PowerVec));
lparam = length(p.Results.params);
for i=1:length(PowerVec)
    FittedParamCell{1,i} = ['c' num2str(PowerVec(i))];
    FittedParamCell{2,i} = PolyCoef(i);
    if i<=lparam
        FittedParamCell(3,i) = p.Results.params(i);
    else
        FittedParamCell(3,i) = FittedParamCell(1,i);
    end
    eval([FittedParamCell{3,i} '=' 'FittedParamCell{2,' num2str(i) '};']);
end

StrFun = [];
for i=1:length(PowerVec)
    if PowerVec(i)>1
        StrFun = [StrFun FittedParamCell{3,i} '*x^' num2str(PowerVec(i)) ''];
    elseif PowerVec(i)==1
        StrFun = [StrFun FittedParamCell{3,i} '*x'];
    else
        StrFun = [StrFun FittedParamCell{3,i}];
    end
    if i<length(PowerVec)
        StrFun = [StrFun ' +'];
    end
end

FittedParamCell(:,end+1) = {'R2'; R2; 'R^{2}'};

span = [min(x) max(x)];
eval(['fun=@(x)' StrFun ';']);
end