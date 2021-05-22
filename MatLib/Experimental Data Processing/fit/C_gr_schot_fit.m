function [str, fun, span] = C_gr_schot_fit(x,y,boxin,f,model,fit_param,varargin)
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

syms vg
independent = {'vg'};
coefficients = {};
coefficients_lim = [];
modelf = f*model;




for p=fit_param
        coefficients{end+1} = char(p{:});
        if ~isempty(fit_cond_val(strcmp(char(p{:}),fit_cond_names),:))
            coefficients_lim(end+1,:) = fit_cond_val(strcmp(char(p{:}),fit_cond_names),:);
        else
        	coefficients_lim(end+1,:) = [-inf inf 1];
        end
end

fo = fitoptions('Method','NonlinearLeastSquares','Lower',coefficients_lim(:,1),'Upper',coefficients_lim(:,2),'StartPoint',coefficients_lim(:,3),'Robust','Bisquare');
str=['Fit Range= '  num2str(boxin)  newline];
ft = fittype(char(modelf), 'independent',independent, 'coefficients',coefficients, 'options',fo);
dvec = ~excludedata(x,y,'box',box);
[fitt, gof] = fit(x(dvec),y(dvec)*f,ft);
%fitt = fit(x,y,ft);
for p=fit_param
    cp = char(p{:});
    str = [str cp '=' num2str(eval(['fitt.' cp]),2) ' '];
    eval([cp '=' num2str(eval(['fitt.' cp])) ';']);
end

str = [str ' R^{2}=' num2str(gof.rsquare,2) ];
%fun = @(vg) A.*sqrt(q.*es.*e0.*N./(2.*n.*(n.*(vb-0.0257)-vg)));
fun = @(vg) 0;
span = [min(x), min(max(x), vb-6*0.0257)];
end

