function [str, fun, span, FittedParamCell] = Ct_Nt_tau_fit(x,y,BoxIn,C0,N,NT,FitProp,varargin)
if isempty(FitProp)
    FitProp = {};
end
dvec = ~isnan(y);
x = x(dvec);
y = y(dvec);
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
%A=1;eox=1;tox=1;es=1;
syms t
independent = {'t'};
coefficients = {};
coefficients_lim = [];
y_fit = y-C0;
dvec = ~excludedata(x,y_fit,'box',Box);
y_fitL = length(y_fit);
y_fitL10 = round(y_fitL/10);
decaydown = mean(y_fit(y_fitL10:2*y_fitL10))>mean(y_fit(end-y_fitL10:end));
y_fit =y_fit/(C0/2);
f = 10^ceil(abs(log10(1/mean(y_fit))));
xd = x(dvec);
xSpan = max(xd) - min(xd);
Cind = 1;
C_t = 0;
if decaydown
    for i=1:NT
        eval(['syms Nt' num2str(i)])
        eval(['syms tau' num2str(i)])
        eval(['C_t = C_t+Nt' num2str(i) '*exp(-t/tau' num2str(i) ');'])
        coefficients{end+1} = ['Nt' num2str(i)];
        coefficients{end+1} = ['tau' num2str(i)];
        cmpN = strcmp(['Nt' num2str(i)],fit_cond_names);
        if any(cmpN)
            coefficients_lim(end+1,:) = fit_cond_val(cmpN,:)/N;
            fit_cond_val(cmpN,:) = [];
            fit_cond_names(cmpN) = [];
        else
            coefficients_lim(end+1,:) = [0 abs(y_fit(Cind)) 10*abs(y_fit(Cind))];
        end
        cmptau = strcmp(['tau' num2str(i)],fit_cond_names);
        if any(cmptau)
            coefficients_lim(end+1,:) = fit_cond_val(cmptau,:);
            fit_cond_val(cmptau,:) = [];
            fit_cond_names(cmptau) = [];
        else
            coefficients_lim(end+1,:) = [0 xSpan*(0.1)^(i-1) Inf];
        end
    end
else
    for i=1:NT
        eval(['syms Nt' num2str(i)])
        eval(['syms tau' num2str(i)])
        eval(['C_t = C_t-Nt' num2str(i) '*exp(-t/tau' num2str(i) ');'])
        coefficients{end+1} = ['Nt' num2str(i)];
        coefficients{end+1} = ['tau' num2str(i)];
        cmpN = strcmp(['Nt' num2str(i)],fit_cond_names);
        if any(cmpN)
            coefficients_lim(end+1,:) = fit_cond_val(cmpN,:)/N;
            fit_cond_val(cmpN,:) = [];
            fit_cond_names(cmpN) = [];
        else
            coefficients_lim(end+1,:) = [0 abs(y_fit(Cind)) 10*abs(y_fit(Cind))];
        end
        cmptau = strcmp(['tau' num2str(i)],fit_cond_names);
        if any(cmptau)
            coefficients_lim(end+1,:) = fit_cond_val(cmptau,:);
            fit_cond_val(cmptau,:) = [];
            fit_cond_names(cmptau) = [];
        else
            coefficients_lim(end+1,:) = [0 xSpan*(0.1)^(i-1) Inf];
        end
    end
end

fo = fitoptions('Method','NonlinearLeastSquares', 'Lower',coefficients_lim(:,1), 'Upper',coefficients_lim(:,3), 'StartPoint',coefficients_lim(:,2), 'TolFun',(1e-16)/f, 'TolX',(1e-16)/f, FitProp{:});
ft = fittype(char(C_t), 'independent',independent, 'coefficients',coefficients, 'options',fo);
[fitt, gof] = fit(x(dvec),y_fit(dvec),ft);

boxstr = [num2str(max(min(x),BoxIn(1)),1) ' to ' num2str(min(max(x),BoxIn(2)),1)];
if length(BoxIn)==4
    boxstr = ['x- ' boxstr ' y- ' num2str(max(min(y),BoxIn(3)),1) ' to ' num2str(min(max(y),BoxIn(4)),1)];
end
str=['Fit Range= ' boxstr  newline];
funstr = '';
FittedParamCell = {};
for i=1:NT
    eval(['Nt' num2str(i) '=' 'N*fitt.Nt' num2str(i) ';']);
    eval(['tau' num2str(i) '=' 'fitt.tau' num2str(i) ';']);
    eval(['Nt=Nt' num2str(i) ';'])
    eval(['tau=tau' num2str(i) ';'])
    str = [str 'Nt_{' num2str(i) '}=' num2str(Nt,2) ' '];
    str = [str '\tau_{' num2str(i) '}=' num2str(tau,2) ' '];
    if decaydown
        funstr = [funstr '+Nt' num2str(i) '*exp(-t/tau' num2str(i) ')'];
    else
        funstr = [funstr '-Nt' num2str(i) '*exp(-t/tau' num2str(i) ')'];
    end
    FittedParamCell(:, end+1) = {['Nt' num2str(i)]; eval(['Nt' num2str(i)]); ['Nt_{' num2str(i) '} - Trap Concentration']};
    FittedParamCell(:, end+1) = {['tau' num2str(i)]; eval(['tau' num2str(i)]); ['\tau_{' num2str(i) '} - Emission Time']};
end

str = [str ' R^{2}=' num2str(gof.rsquare,3) ];
eval(['fun = @(t) (' funstr ')*(C0/2/N)+C0;']);
span = [min(x), max(x)];
end

