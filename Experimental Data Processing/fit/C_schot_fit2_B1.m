function [str, fun, span, FittedVarsCell] = C_schot_fit2_B1(x,y,A,es,N,Vb,n,a,n0,varargin)

p = inputParser;
p.KeepUnmatched=true;
p.addParameter('range', [-inf inf -inf inf], @isnumeric);
p.addParameter('FirstOut', '', @isstring);
p.addParameter('T', 296, @isnumeric);
p.addParameter('f', 10^ceil(abs(log10(1/mean(y)))), @isnumeric);
p.addParameter('N',[0 1e16 inf], @isnumeric);
p.addParameter('Vb',[0 0.5 10], @isnumeric);
p.addParameter('n',[0.01 1 1000], @isnumeric);
p.addParameter('n0',[0 9e12 inf], @isnumeric);
p.addParameter('a',[-10 0.5 10 ], @isnumeric);
p.parse(varargin{:});

T = p.Results.T;
f = p.Results.f;

if length(p.Results.range)==2
    range = [p.Results.range -inf inf];
else
    range = p.Results.range;
end

syms V
independent = {'V'};
k = 8.617e-5; % [eV/T]
h = 4.1356e-15;
q = 1.6e-19; % [C]
e0 = 8.854e-14;
sqef = sqrt(q*e0*es)*f;
v_fermi_gr = 1e8; % [cm/sec]
coefficients = {};
coefficients_lim = [];

if isa(N,'sym')
    syms Nsq
    coefficients{end+1} = 'Nsq';
    coefficients_lim(end+1,:) = sqrt(p.Results.N);
else 
    Nsq = sqrt(N);
end

if isa(Vb,'sym')
    coefficients{end+1} = 'Vb';
    coefficients_lim(end+1,:) = p.Results.Vb;
end
if isa(n,'sym')
    coefficients{end+1} = 'n';
    coefficients_lim(end+1,:) = p.Results.n;
end

if isa(n0,'sym')
    syms n0sq
    coefficients{end+1} = 'n0sq';
    coefficients_lim(end+1,:) = sqrt(p.Results.n0);
else 
    n0sq = sqrt(n0);
end

if isa(a,'sym')
    coefficients{end+1} = 'a';
    coefficients_lim(end+1,:) = p.Results.a;
else
    a = h/4/sqrt(pi)*v_fermi_gr*Nsq*sqrt(es*e0/2/q)/n0sq;
end


cschot = A*sqef*Nsq/sqrt(2*n*(n*(Vb-k*T)-V+a*V));


fo = fitoptions('Method','NonlinearLeastSquares','Lower',coefficients_lim(:,1),'Upper',coefficients_lim(:,3),'StartPoint',coefficients_lim(:,2),'Robust','Bisquare');
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

