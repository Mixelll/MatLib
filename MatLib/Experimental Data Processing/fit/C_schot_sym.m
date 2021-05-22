function [cschot, varargout] = C_schot_sym(A,es,N,vb,varargin)
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
cschot = A/(sqrt(2*(vg-vb)/qes/N));
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

