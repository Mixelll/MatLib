function [str, fun, span, FittedVarsCell] = C_schot_fit2_B1(x,y,A,es,N,Vb,n,a,n0,varargin)
p = inputParser;
p.KeepUnmatched=true;
p.addParameter('Range', [-inf inf], @isnumeric);
p.addParameter('T', 296, @isnumeric);
p.addParameter('f', 10^ceil(abs(log10(1/mean(y)))), @isnumeric);
p.addParameter('Doping_Lim',[0 1e17 inf], @isnumeric);
p.addParameter('Vb_Lim',[0 0.5 10], @isnumeric);
p.addParameter('Ideality_Factor_Lim',[0.01 1 1000], @isnumeric);
p.addParameter('Initial_Graphene_Doping_Lim',[0 9e12 inf], @isnumeric);
p.addParameter('v_fermi_Gr',1e8, @isnumeric);
p.addParameter('Graphene_Vdrop_Coeff_Lim',[-10 0.5 10 ], @isnumeric);
p.addParameter('FitProperties',{'Robust','Bisquare'});
p.addParameter('FirstOut', '', @isstring);
p.parse(varargin{:});
if strcmpi(x, 'parser')
    s = struct;
    for c = p.UsingDefaults
        s.(c{:}) = p.Results.(c{:});
    end
    str = s; fun = []; span = []; FittedVarsCell = [];
elseif strcmpi(x, 'model')
    str = 'C = A*sqrt(es*e0*q*N)/sqrt(2*n*(n*(Vb-k*T)-V+a*V)), a = h/4/sqrt(pi)*v_fermi_gr*sqrt(es*e0*N)/sqrt(2*q*n0)'; fun = []; span = []; FittedVarsCell = [];
else
T = p.Results.T;
f = p.Results.f;
FitProp = p.Results.FitProperties;
if isempty(FitProp)
    FitProp = {};
elseif ~iscell(FitProp)
    FitProp = {FitProp};
end

if length(p.Results.Range)==2
    range = [p.Results.Range -inf inf];
else
    range = p.Results.Range;
end

syms V
independent = {'V'};
k = 8.617e-5; % [eV/T]
h = 4.1356e-15;
q = 1.6e-19; % [C]
e0 = 8.854e-14;
sqef = sqrt(q*e0*es)*f;
v_fermi_gr = p.Results.v_fermi_Gr; % [cm/sec]
coefficients = {};
coefficients_lim = [];

if isa(N,'sym')
    syms Nsq
    coefficients{end+1} = 'Nsq';
    coefficients_lim(end+1,:) = sqrt(p.Results.Doping_Lim);
else 
    Nsq = sqrt(N);
end

if isa(Vb,'sym')
    coefficients{end+1} = 'Vb';
    coefficients_lim(end+1,:) = p.Results.Vb_Lim;
end
if isa(n,'sym')
    coefficients{end+1} = 'n';
    coefficients_lim(end+1,:) = p.Results.Ideality_Factor_Lim;
end

if isa(n0,'sym')
    syms n0sq
    coefficients{end+1} = 'n0sq';
    coefficients_lim(end+1,:) = sqrt(p.Results.Initial_Graphene_Doping_Lim);
else 
    n0sq = sqrt(n0);
end

if isa(a,'sym')
    coefficients{end+1} = 'a';
    coefficients_lim(end+1,:) = p.Results.Graphene_Vdrop_Coeff_Lim;
else
    a = h/4/sqrt(pi)*v_fermi_gr*Nsq*sqrt(es*e0/2/q)/n0sq;
end


cschot = A*sqef*Nsq/sqrt(2*n*(n*(Vb-k*T)-V+a*V));


fo = fitoptions('Method','NonlinearLeastSquares','Lower',coefficients_lim(:,1),'Upper',coefficients_lim(:,3),'StartPoint',coefficients_lim(:,2), FitProp{:});
ft = fittype(char(cschot), 'independent',independent, 'coefficients',coefficients, 'options',fo);
dvec = ~excludedata(x,y,'box',range) & ~isnan(y);
Vbias = x(dvec);
C = y(dvec);
[fitt, gof] = fit(Vbias,C*f,ft);

if abs(range(1)-Vbias(1))<0.1, x1 = range(1); else, x1 = Vbias(1); end
if abs(range(2)-Vbias(end))<0.1, x2 = range(2); else, x2 = Vbias(end); end
box = [x1 x2];
str=['Fit Range= '  num2str(box,2)  newline];
FittedVarsCell = {};
if any(strcmp(coefficients, 'Nsq'))
    eval(['N=' num2str(fitt.Nsq^2) ';']);
    FittedVarsCell(:, end+1) = {'N'; N; 'N - Carrier Density [cm^{-3}]'};
end
str = [str 'N=' num2str(N,2) ' '];
if any(strcmp(coefficients, 'Vb'))
    eval(['Vb=' num2str(fitt.Vb) ';']);
    FittedVarsCell(:, end+1) = {'Vb'; Vb; 'Vb - V Built in [V]'};
end
str = [str 'Vb=' num2str(Vb,2)];
if any(strcmp(coefficients, 'n'))
    eval(['n=' num2str(fitt.n) ';']);
    FittedVarsCell(:, end+1) = {'n'; n; 'n - Ideality Factor'};
end
str = [str ' n=' num2str(n,2)];
if any(strcmp(coefficients, 'n0sq'))
    eval(['n0=' num2str(fitt.n0sq^2) ';']);
    FittedVarsCell(:, end+1) = {'n0'; n0; 'n0 - Graphene Doping [cm^{-2}]'};
end
str = [str ' n0=' num2str(n0,2)];
if any(strcmp(coefficients, 'a'))
    eval(['a=' num2str(fitt.a) ';']);
else
    a = h/4/sqrt(pi)*v_fermi_gr*sqrt(es*e0*N/2/q/n0);
end
FittedVarsCell(:, end+1) = {'a'; a; 'a - Graphene Ef Shift Param'};
str = [str ' a=' num2str(a,2)];
str = [str ' R^{2}=' num2str(gof.rsquare,3) ];
FittedVarsCell(:, end+1) = {'R2'; gof.rsquare; 'R^{2}'};
fun = @(V) A.*sqrt(q.*es.*e0.*N./(2.*n.*(n.*(Vb-k*T)-V +a.*V)));
span = [min(x), min(max(x), Vb-6*k*T)];
if ~isempty(p.ResultsFirstOut)
    str = cellfun(@(c) FittedVarsCell(2,strcmpi(c, FittedVarsCell(1,:))),p.ResultsFirstOut);
end
end
end

