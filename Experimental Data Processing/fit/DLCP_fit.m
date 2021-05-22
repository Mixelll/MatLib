function [str, fun, span, FittedParamCell] = DLCP_fit(x,y, A, varargin)
p = inputParser;
p.addParameter('range', [-inf inf], @isnumeric);
p.addParameter('FirstOut', '', @isstring);
p.addParameter('T', 296);
% p.addParameter('f', 10^ceil(abs(log10(1/y(end-1)))), @isnumeric);
% p.addParameter('N',[0 1e16 inf], @isnumeric);
p.addParameter('es',11.68);
p.addParameter('Order',1);
p.parse(varargin{:});

% T = p.Results.T;
Order = p.Results.Order;
if ~isnumeric(Order)
    Order = 1;
end
if length(p.Results.range)==2
    Range = [p.Results.range -inf inf];
else
    Range = p.Results.range;
end

% k = 8.617e-5;
% h = 4.1356e-15;
q = 1.6e-19;
e0 = 8.854e-14;
es = p.Results.es;
if ~isnumeric(es)
    es = 11.68;
end
% v_fermi_gr = 1e8; %cm/sec

params = {'C0','C1'};
if Order>1
    params = {'C0' 'C1' 'C2'};
end
[~, fun, span, FittedParamCell] = str_poly_fit(x,y,Order, 'params',params, 'range',Range);
C0 = FittedParamCell{2, strcmpi(FittedParamCell(3,:),params{1})};
C1 = FittedParamCell{2, strcmpi(FittedParamCell(3,:),params{2})};
if Order>1
    C1 = FittedParamCell{2, strcmpi(FittedParamCell(3,:),params{3})};
end
N = -C0^3 / (2*q*e0*es*A^2*C1);

str='';
for c = FittedParamCell
    str = [str ', ' c{3} '=' num2str(c{2},2)];
end
str(1:2) = '';
FittedParamCell = [{'N'; N;  'N - Carrier Density [cm^{-3}]'} FittedParamCell];
str = ['N - Carrier Density = ' num2str(FittedParamCell{2,1},2) 'cm^{-3}' newline str];

Boxstr = [num2str(max(min(x(1),x(end)),Range(1)),2) ' to ' num2str(min(max(x(1),x(end)),Range(2)),2)];
if length(p.Results.range)==4
    Boxstr = ['x- ' Boxstr ' y- ' num2str(max(min(y),Range(3)),2) ' to ' num2str(min(max(y),Range(4)),2)];
end
str=['Fit Range= '  num2str(Boxstr)  newline str];
end

