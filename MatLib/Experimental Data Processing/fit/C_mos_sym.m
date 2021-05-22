function [cmos, varargout] = C_mos_sym(A,eox,tox,es,N,vfb,varargin)
if ~isempty(varargin)
    if isa(varargin{1}, 'numeric')
        fit_cond_val = cell2mat(varargin(1:2:end)');
        fit_cond_names = varargin(2:2:end);
    else
        fit_cond_val = cell2mat(varargin(2:2:end)');
        fit_cond_names = varargin(1:2:end);
    end
end

%A=1;eox=1;tox=1;es=1;
syms vg
independent = {'vg'};
q = 1.6e-19;
e0 = 8.854e-14;
qes = q*es*e0;
cox = eox*e0/tox;
v1 = 2*(vg -vfb);
v2 = 2*qes*N/cox^2;
fis = (v1 + v2 + sqrt(v2^2 + 2*v1*v2))/2;
cmos = A/(1/cox + sqrt(2*fis/qes/N));
coefficients = {};
problem = {};
for i={A,eox,tox,es,N,vfb}
    if isa(i,'sym')
        coefficients{end+1} = char(i);
        coefficients_lim(end+1,:) = fit_cond_val(strcmp(char(i),fit_cond_names),:);
    else
        problem{end+1} = char(i);
    end
end   
varargout = {independent, coefficients, problem, coefficients_lim};
end

