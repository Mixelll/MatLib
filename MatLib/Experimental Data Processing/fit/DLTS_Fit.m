function [str, fun, span, FittedParamCell] = DLTS_Fit(T,e1,varargin)
p = inputParser;
p.KeepUnmatched=true;
p.addParameter('Range', [-inf inf], @isnumeric); % Use solve to extract tau
p.parse(varargin{:});
Range = p.Results.Range;

kB = 8.617e-5;

if size(Range,2)==2
    Box = [Range -inf inf];
elseif size(Range,1)==2
    Box = [Range.' -inf inf];
elseif isempty(Range)
    Box = [-inf inf -inf inf];
end

dvec = isnan(e1);
if any(dvec)
    T = T(~dvec);
    e1 = e1(~dvec);
end

dvec = ~excludedata(T,e1,'box',Box);
T1 = 1000./T(dvec);
e1T2 = log(e1(dvec)./T(dvec).^2);
[OF, S] = polyfit(T1, e1T2, 1); E = -1000*kB*OF(1); A = exp(OF(2)); R2 = 1 - (S.normr/norm(e1T2 - mean(e1T2)))^2;

boxstr = [num2str(max(min(T),Range(1)),3) ' to ' num2str(min(max(T),Range(2)),3)];
if length(Range)==4
    boxstr = [boxstr ' e1- ' num2str(max(min(e1),Range(3)),3) ' to ' num2str(min(max(e1),Range(4)),3)];
end
str=['Fit Range T[K] = ' boxstr  newline];

str = [str 'E = ' num2str(E,3) ' [eV], A = ' num2str(A,3) ' [sec^{-1}], R^{2} = ' num2str(R2,3)];
fun = @(x) A*exp(-E*x/(1000*kB));
span = [min(T1), max(T1)];
FittedParamCell = {'E' ; E ;'Trap Activation Energy [eV]'};
FittedParamCell(:, end+1) = {'A' ; A ; 'Temperature Independent pre-exponential Factor [sec^{-1}]'};
FittedParamCell(:, end+1) = {'R' ; R2 ; 'R^{2}'};

end

