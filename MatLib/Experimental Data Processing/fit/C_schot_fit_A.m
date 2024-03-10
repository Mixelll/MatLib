function [str, fun, span, FittedParamCell] = C_schot_fit_A(x, y, BoxIn, A, es, N, Vb, n, T, FitProp, varargin)
if strcmpi(x, 'model')
    str = 'C = A*sqrt(es*e0*q*N)/sqrt(2*n*(n*(Vb-k*T)-V))'; fun = []; span = []; FittedParamCell = [];
else
NotNumSym = @(x) ~isnumeric(x) && ~isa(x,'sym');
if isempty(FitProp)
    FitProp = {};
elseif ~iscell(FitProp)
    FitProp = {FitProp};
end
if size(x,2)>size(x,1)
    x = x.';
end
if size(y,2)>size(y,1)
    y = y.';
end
dvec = ~isnan(y) & ~isnan(x);
x = x(dvec);
y = y(dvec);
f = 10^ceil(abs(log10(1/mean(y))));
if ~isempty(varargin)
    if isa(varargin{1}, 'numeric')
        fit_cond_val = cell2mat(varargin(1:2:end)');
        fit_cond_names = varargin(2:2:end);
    else
        fit_cond_val = cell2mat(varargin(2:2:end)');
        fit_cond_names = varargin(1:2:end);
    end
else
    fit_cond_val = [];
    fit_cond_names = {};
end
if size(BoxIn,2)==2
    Box = [BoxIn -inf inf];
elseif size(BoxIn,1)==2
    Box = [BoxIn.' -inf inf];
elseif isempty(BoxIn)
    Box = [-inf inf -inf inf];
end

syms c V
independent = {'V'};
if ~T
    T = 295;
end
k = 8.617e-5;
q = 1.6e-19;
e0 = 8.854e-14;
sqef = sqrt(q*e0)*f;
cschot = c/sqrt(2*n*(n*(Vb-0.0257)-V));
coefficients = {};
sym_factors = {};
sym_titles = {};
sym_title_find = {'N'; 'N - Carrier Density [cm^{-3}]'};
coefficients_lim = [];
prod_values = [];


if NotNumSym(es), es = sym('es'); end
if NotNumSym(N), N = sym('N'); end
if NotNumSym(Vb), Vb = sym('Vb'); end
if NotNumSym(n), n = sym('n'); end

if isa(Vb,'sym')
    coefficients{end+1} = 'Vb';
    cmp = strcmp('Vb',fit_cond_names);
    if any(cmp)
        coefficients_lim(end+1,:) = fit_cond_val(cmp,:);
        fit_cond_val(cmp,:) = [];
        fit_cond_names(cmp) = [];
    else
        coefficients_lim(end+1,:) = [0 0.5 10];
    end
end
if isa(n,'sym')
    coefficients{end+1} = 'n';
    cmp = strcmp('n',fit_cond_names);
    if any(cmp)
        coefficients_lim(end+1,:) = fit_cond_val(cmp,:);
        fit_cond_val(cmp,:) = [];
        fit_cond_names(cmp) = [];
    else
        coefficients_lim(end+1,:) = [0.01 1 100];
    end
end



sqrt_param = {'es','N'};
char_param = {'A','es','N'};

% only linear param 
param = {A,es,N};
flag = 1;
for i=1:length(param)
    p = param{i};
    if isa(p,'sym')
        if flag==1
            coefficients{end+1} = 'c';
            coefficients_lim(end+1,:) = [sqef sqef sqef];
            flag = 0;
        end
        sym_factors{end+1} = char(p);
        sym_titles(end+1) = sym_title_find(2,strcmpi(sym_title_find(1,:),char(p)));
        if ~isempty(fit_cond_val(strcmp(char(p),fit_cond_names),:))
            if any(strcmp(char(p),sqrt_param))
                coefficients_lim(end,:) = coefficients_lim(end,:).*sqrt(abs(fit_cond_val(strcmp(char(p),fit_cond_names),:)));
            else
                coefficients_lim(end,:) = coefficients_lim(end,:).*fit_cond_val(strcmp(char(p),fit_cond_names),:);
            end
        else
        	coefficients_lim(end,:) = [0 1 inf];
        end
    else
        if any(strcmp(char_param{i},sqrt_param))
            if flag==1
                coefficients{end+1} = 'c';
                coefficients_lim(end+1,:) = [sqef sqef sqef]*sqrt(p);
                flag = 0;
            else
                coefficients_lim(end,:) = coefficients_lim(end,:)*sqrt(p);
            end
            prod_values(end+1) = sqrt(p);
        else
            if flag==1
                coefficients{end+1} = 'c';
                coefficients_lim(end+1,:) = [sqef sqef sqef]*p;
                flag = 0;
            else
                coefficients_lim(end,:) = coefficients_lim(end,:)*p;
            end
            prod_values(end+1) = p;
        end
    end
end
% only linear param ^

Boxstr = [num2str(max(min(x(1),x(end)),BoxIn(1)),2) ' to ' num2str(min(max(x(1),x(end)),BoxIn(2)),2)];
if length(BoxIn)==4
    Boxstr = ['x- ' Boxstr ' y- ' num2str(max(min(y),BoxIn(3)),2) ' to ' num2str(min(max(y),BoxIn(4)),2)];
end
fo = fitoptions('Method','NonlinearLeastSquares', 'Lower',coefficients_lim(:,1), 'StartPoint',coefficients_lim(:,2), 'Upper',coefficients_lim(:,3), FitProp{:});
str=['Fit Range= '  num2str(Boxstr)  newline];
if coefficients{end} ~= 'c'
    c = A*sqrt(N)*sqef;
end
ft = fittype(char(cschot), 'independent',independent, 'coefficients',coefficients, 'options',fo);
dvec = ~excludedata(x,y,'box',Box);
try
    [fitt, gof] = fit(x(dvec),y(dvec)*f,ft);
catch ME
    if contains(ME.identifier, 'complex','IgnoreCase',true)
        warning('Complex value computed by model function, fitting cannot continue. Assisgning NaN values to outputs.');
        fitt = [];
        gof = [];
    end
end
        
if ~isempty(fitt)
    FittedParamCell = {};
    if coefficients{end} == 'c' && ~isempty(sym_factors)
        if any(strcmp(sym_factors{1},sqrt_param))
            str = [str sym_factors{:} '=' num2str((fitt.c/sqef/prod(prod_values))^2,2) ' '];
            eval([sym_factors{1} '=' num2str((fitt.c/sqef/prod(prod_values))^2) ';']);
            FittedParamCell(:, end+1) = {sym_factors{1}; (fitt.c/sqef/prod(prod_values))^2; sym_titles{1}};
        else
            str = [str sym_factors{:} '=' num2str(fitt.c/sqef/prod(prod_values),2) ' '];
            eval([sym_factors{1} '=' num2str(fitt.c/sqef/prod(prod_values)) ';']);
            FittedParamCell(:, end+1) = {sym_factors{1}; fitt.c/sqef/prod(prod_values); sym_titles{1}};
        end
    end
    if any(strcmp(coefficients, 'Vb'))
        str = [str 'Vb=' num2str(fitt.Vb,2)];
        eval(['Vb=' num2str(fitt.Vb) ';']);
        FittedParamCell(:, end+1) = {'Vb'; Vb; 'Vb - V Built in [V]'};
    end
    if any(strcmp(coefficients, 'n'))
        str = [str ' n=' num2str(fitt.n,2)];
        eval(['n=' num2str(fitt.n) ';']);
        FittedParamCell(:, end+1) = {'n'; n; 'n - Ideality Factor'};
    end
    str = [str ' R^{2}=' num2str(gof.rsquare,3) ];
    FittedParamCell(:, end+1) = {'R2'; gof.rsquare; 'R^{2}'};
    fun = @(V) A.*sqrt(q.*es.*e0.*N./(2.*n.*(n.*(Vb-k*T)-V)));
    span = [min(x), min(max(x), Vb-2*k*T)];
else
    if coefficients{end} == 'c' && ~isempty(sym_factors)
        coefficients(end)= [];
        coefficients = [sym_factors(1) coefficients];
    end
    str = ['Complex Values Computed.' newline 'Fitting aborted'];
    fun = [];
    span = [];
    FittedParamCell = [coefficients {'R2'}; num2cell(nan(1,length(coefficients)+1)); coefficients {'R2'}]; 
end
end

