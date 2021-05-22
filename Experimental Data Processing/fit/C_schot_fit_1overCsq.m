function [str, fun, span] = C_schot_fit_1overCsq(x,y,boxin,f,over,A,es,N,Vb,n,varargin)
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
if length(boxin)==2
    box = [boxin -inf inf];
elseif isempty(boxin)
    box = [-inf inf -inf inf];
end
%A=1;eox=1;tox=1;es=1;
syms c V
independent = {'V'};
T = 295;
k = 8.617e-5;
q = 1.6e-19;
e0 = 8.854e-14;
qef = q*e0/f;
cschot = (2*n*(n*(Vb-0.0257)-V))/c;
coefficients = {};
sym_factors = {};
coefficients_lim = [];
prod_values = [];


if isa(Vb,'sym')
    coefficients{end+1} = 'Vb';
    cmp = strcmp('Vb',fit_cond_names);
    if any(cmp)
        coefficients_lim(end+1,:) = fit_cond_val(cmp,:);
        fit_cond_val(cmp,:) = [];
        fit_cond_names(cmp) = [];
    else
        coefficients_lim(end+1,:) = [-inf inf 0];
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
        coefficients_lim(end+1,:) = [-inf inf 0];
    end
end



p2_param = {'A'};
char_param = {'A','es','N'};

% only linear param 
param = {A,es,N};
flag = 1;
for i=1:length(param)
    p = param{i};
    if isa(p,'sym')
        if flag==1
            coefficients{end+1} = 'c';
            coefficients_lim(end+1,:) = [qef qef qef];
            flag = 0;
        end
        sym_factors{end+1} = char(p);
        if ~isempty(fit_cond_val(strcmp(char(p),fit_cond_names),:))
            if any(strcmp(char(p),p2_param))
                coefficients_lim(end,:) = coefficients_lim(end,:).*abs(fit_cond_val(strcmp(char(p),fit_cond_names),:)).^2;
            else
                coefficients_lim(end,:) = coefficients_lim(end,:).*fit_cond_val(strcmp(char(p),fit_cond_names),:);
            end
        else
        	coefficients_lim(end,:) = [-inf inf 0];
        end
    else
        if any(strcmp(char_param{i},p2_param))
            if flag==1
                coefficients{end+1} = 'c';
                coefficients_lim(end+1,:) = [qef qef qef]*p.^2;
                flag = 0;
            else
                coefficients_lim(end,:) = coefficients_lim(end,:)*p.^2;
            end
            prod_values(end+1) = p.^2;
        else
            if flag==1
                coefficients{end+1} = 'c';
                coefficients_lim(end+1,:) = [qef qef qef]*p;
                flag = 0;
            else
                coefficients_lim(end,:) = coefficients_lim(end,:)*p;
            end
            prod_values(end+1) = p;
        end
    end
end
% only linear param ^

fo = fitoptions('Method','NonlinearLeastSquares','Lower',coefficients_lim(:,1),'Upper',coefficients_lim(:,2),'StartPoint',coefficients_lim(:,3),'Robust','Bisquare');
str=['Fit Range= '  num2str(boxin)  newline];
if coefficients{end} ~= 'c'
    c = A*sqrt(N)*qef;
end
ft = fittype(char(cschot), 'independent',independent, 'coefficients',coefficients, 'options',fo);
dvec = ~excludedata(x,y,'box',box);
y = 1./y.^2;
[fitt, gof] = fit(x(dvec),y(dvec)*f,ft);
%fitt = fit(x,y,ft);
if coefficients{end} == 'c' && ~isempty(sym_factors)
    if any(strcmp(sym_factors{1},p2_param))
        str = [str sym_factors{:} '=' num2str((fitt.c/qef/sqrt(prod(prod_values))),2) ' '];
        eval([sym_factors{1} '=' num2str((fitt.c/qef/sqrt(prod(prod_values)))) ';']);
    else
        str = [str sym_factors{:} '=' num2str(fitt.c/qef/prod(prod_values),2) ' '];
        eval([sym_factors{1} '=' num2str(fitt.c/qef/prod(prod_values)) ';']);
    end
end
if any(strcmp(coefficients, 'Vb'))
    str = [str 'Vb=' num2str(fitt.Vb,2)];
    eval(['Vb=' num2str(fitt.Vb) ';']);
end
if any(strcmp(coefficients, 'n'))
    str = [str ' n=' num2str(fitt.n,2)];
    eval(['n=' num2str(fitt.n) ';']);
end
str = [str ' R^{2}=' num2str(gof.rsquare,3) ];
if ~over
    fun = @(V) A.*sqrt(q.*es.*e0.*N./(2.*n.*(n.*(Vb-k*T)-V)));
else
    fun = @(V) (2.*n.*(n.*(Vb-k*T)-V))./(q.*es.*e0.*N.*A.^2);
end
span = [min(x), min(max(x), Vb-6*k*T)];
end

