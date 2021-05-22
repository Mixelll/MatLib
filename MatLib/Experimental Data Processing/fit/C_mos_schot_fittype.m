function ft = C_mos_schot_fittype(Ac,eox,tox,es,N,vfb,As,vb,varargin)
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

%A=1;eox=1;tox=1;es=1;
syms vg
independent = {'vg'};
q = 1.6e-19;
e0 = 8.854e-14;
qes = q*es*e0;
cox = eox*e0/tox;
v1 = 2*abs(vg-vfb);
v2 = 2*qes*N/cox^2;
fis = (v1 + v2 + sqrt(v2^2 + 2*v1*v2))/2;
ctot = Ac/(1/cox + sqrt(2*fis/qes/N)) + As/(sqrt(2*abs(vg-vb)/qes/N));
disp(ctot)
coefficients = {};
coefficients_lim = [];
problem = {};

char_param = {'Ac','eox','tox','es','N','vfb','As','vb', 'n'};
param = {Ac,eox,tox,es,N,vfb,As,vb};
for i=1:length(param)
    p = param{i};
    if isa(p,'sym')&&contains(char(ctot),char(p))
        coefficients{end+1} = char(p);
        if ~isempty(fit_cond_val(strcmp(char(p),fit_cond_names),:))
            coefficients_lim(end+1,:) = fit_cond_val(strcmp(char(p),fit_cond_names),:);
        else
            coefficients_lim(end+1,:) = [-inf inf 0];
        end
    else
        problem(end+1) = char_param(i);
    end
end
if ~isempty(coefficients_lim)
    fo = fitoptions('Method','NonlinearLeastSquares','Lower',coefficients_lim(:,1),'Upper',coefficients_lim(:,2),'StartPoint',coefficients_lim(:,3),'DiffMaxChange',1e80,'DiffMinChange',1e-100,'TolFun',1e-100,'TolX',1e-100,'MaxFunEvals',10000,'MaxIter',10000,'robust','LAR');
else
    fo = fitoptions('Method','NonlinearLeastSquares','DiffMaxChange',1e80,'DiffMinChange',1e-80,'TolFun',1e-80,'TolX',1e-80,'MaxFunEvals',5000,'MaxIter',5000);
end
disp(fo)
disp(coefficients)
ft = fittype(char(ctot), 'independent',independent, 'coefficients',coefficients, 'options',fo);
end

