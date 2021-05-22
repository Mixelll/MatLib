function [str, fun, span, fittedvar] = I_gr_si_schot_forward_Trushin(x,y,boxin,f,T,R,A,VB,n,n0,outfittedvar,varargin)
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
if isempty(T)
    T = 296;
end
k = 8.617e-5;
h = 4.1356e-15;
h_bar = h/(2*pi);
q = 1.6e-19;
e0 = 8.854e-14;
v_fermi = 1e8; %m/sec
coefficients = {};
coefficients_lim = [];


if isa(VB,'sym')
    coefficients{end+1} = 'VB';
    cmp = strcmp('VB',fit_cond_names);
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
    syms n0sq
    coefficients{end+1} = 'n0sq';
    cmp = strcmp('n0',fit_cond_names);
    if any(cmp)
        coefficients_lim(end+1,:) =  sqrt(fit_cond_val(cmp,:));
        fit_cond_val(cmp,:) = [];
        fit_cond_names(cmp) = [];
    else
        coefficients_lim(end+1,:) = [0 3e6 inf];
    end
else 
    n0sq = sqrt(n0);
end

E_F = h/2/sqrt(pi)*v_fermi*n0sq;
Astar = (VB + k*T)*k^3/(pi^2*h_bar^3*v_fermi^2*E_F);
I0 = A*Astar*T^3*exp(-VB/k/T);
Iforward = I0*(exp(V/(n*k*T))-1)*f;

fo = fitoptions('Method','NonlinearLeastSquares','Lower',coefficients_lim(:,1),'Upper',coefficients_lim(:,3),'StartPoint',coefficients_lim(:,2),'Robust','Bisquare');
ft = fittype(char(Iforward), 'independent',independent, 'coefficients',coefficients, 'options',fo);
dvec = ~excludedata(x,y,'box',box);
Vbias = x(dvec);
I = y(dvec);
[fitt, gof] = fit(Vbias-I*R,I*f,ft);
fittedvar = [];

[~, ay] = min(abs(x-box(1)));
[~, by] = min(abs(x-box(2)));
str = ['Fit Range: x =['  num2str([box(1)-R*y(ay) box(2)-R*y(by)]) ']  y=[' num2str(box(3:4)) ']'  newline];

if any(strcmp(coefficients, 'VB'))
    eval(['VB=' num2str(fitt.VB) ';']);
    fittedvar = [fittedvar VB];
end
str = [str 'VB=' num2str(VB,3)];
if any(strcmp(coefficients, 'n'))
    eval(['n=' num2str(fitt.n) ';']);
    fittedvar = [fittedvar n];
end
str = [str ' n=' num2str(n,2)];
if any(strcmp(coefficients, 'n0sq'))
    eval(['n0=' num2str(fitt.n0sq^2) ';']);
    fittedvar = [fittedvar n0];
end
str = [str ' n0=' num2str(n0,2)];

E_F = h/2/sqrt(pi)*v_fermi*sqrt(n0);
str = [str ' Ef=' num2str(E_F,2)];

str = [str ' R^{2}=' num2str(gof.rsquare,3) ];
fun = @(V) A*(VB + k*T)*k^3*T^3./(pi^2*h_bar^3*v_fermi^2*E_F).*exp(-VB/k/T).*(exp(V./(n*k*T))-1);
span = [min(0,box(1)-R*y(ay)) max(max(x-R*y),box(2)-R*y(by))];
if ~isempty(outfittedvar)
    str = fittedvar(outfittedvar);
end
end

