function [str, fun, span] = C_schot_fit2_BG1(x,y,funsqrt,boxin,f,A,es,N,Vb,n,a,n0,varargin)
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
syms V
independent = {'V'};
T = 295;
k = 8.617e-5;
h = 4.1356e-15;
q = 1.6e-19;
e0 = 8.854e-14;
sqef = sqrt(q*e0*es)*f;
v_fermi = 1e8; %m/ses
coefficients = {};
coefficients_lim = [];

if isa(N,'sym')
    syms Nsq
    coefficients{end+1} = 'Nsq';
    cmp = strcmp('N',fit_cond_names);
    if any(cmp)
        coefficients_lim(end+1,:) = sqrt(fit_cond_val(cmp,:));
        fit_cond_val(cmp,:) = [];
        fit_cond_names(cmp) = [];
    else
        coefficients_lim(end+1,:) = [0 0.5 10];
    end
else 
    Nsq = sqrt(N);
end
if isa(Vb,'sym')
    coefficients{end+1} = 'Vb';
    cmp = strcmp('Vb',fit_cond_names);
    if any(cmp)
        coefficients_lim(end+1,:) = fit_cond_val(cmp,:);
        fit_cond_val(cmp,:) = [];
        fit_cond_names(cmp) = [];
    else
        coefficients_lim(end+1,:) = [-10 0.5 10];
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
        coefficients_lim(end+1,:) = [0.01 1 1000];
    end
end
if isa(n0,'sym')
    coefficients{end+1} = 'n0';
    cmp = strcmp('n0',fit_cond_names);
    if any(cmp)
        coefficients_lim(end+1,:) = fit_cond_val(cmp,:);
        fit_cond_val(cmp,:) = [];
        fit_cond_names(cmp) = [];
    else
        coefficients_lim(end+1,:) = [0 1e13 inf];
    end
end
if isa(a,'sym')
    coefficients{end+1} = 'a';
    cmp = strcmp('a',fit_cond_names);
    if any(cmp)
        coefficients_lim(end+1,:) = fit_cond_val(cmp,:);
        fit_cond_val(cmp,:) = [];
        fit_cond_names(cmp) = [];
    else
        coefficients_lim(end+1,:) = [-10 0.5 10 ];
    end
else
    a = h/4/sqrt(pi)*v_fermi*Nsq*sqrt(es*e0/2/q)/sqrt(n0);
end


cschot = A*sqef*Nsq/funsqrt(n,vb,a,V);


fo = fitoptions('Method','NonlinearLeastSquares','Lower',coefficients_lim(:,1),'Upper',coefficients_lim(:,3),'StartPoint',coefficients_lim(:,2),'Robust','Bisquare');
str=['Fit Range= '  num2str(boxin)  newline];
if coefficients{end} ~= 'c'
    c = A*sqrt(N)*sqef;
end
ft = fittype(char(cschot), 'independent',independent, 'coefficients',coefficients, 'options',fo);
dvec = ~excludedata(x,y,'box',box);
[fitt, gof] = fit(x(dvec),y(dvec)*f,ft);
%fitt = fit(x,y,ft);
if any(strcmp(coefficients, 'Nsq'))
    eval(['N=' num2str(fitt.Nsq^2) ';']);
end
str = [str 'N=' num2str(Nsq^2,2) ' '];
if any(strcmp(coefficients, 'Vb'))
    eval(['Vb=' num2str(fitt.Vb) ';']);
end
str = [str 'Vb=' num2str(Vb,2)];
if any(strcmp(coefficients, 'n'))
    eval(['n=' num2str(fitt.n) ';']);
end
str = [str ' n=' num2str(n,2)];
if any(strcmp(coefficients, 'n0'))
    eval(['n0=' num2str(fitt.n0) ';']);
end
str = [str ' n0=' num2str(n0,2)];
if any(strcmp(coefficients, 'a'))
    eval(['a=' num2str(fitt.a) ';']);
else
    a = h/4/sqrt(pi)*v_fermi*sqrt(es*e0*N/2/q/n0);
end
str = [str ' a=' num2str(a,2)];
str = [str ' R^{2}=' num2str(gof.rsquare,3) ];
fun = @(V) A.*sqrt(q.*es.*e0.*N./funsqrt(n,vb,a,V);
span = [min(x), min(max(x), Vb-6*k*T)];
end

