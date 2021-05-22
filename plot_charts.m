function [s_diffdc, s_diffdc_integ] = plot_charts(scells, line_plot, integplot, smooths_init, dcs, ref, plotop1, plotop2, trimleg, fit_data, headerparam, boxin)

if length(boxin)==2
    box = [boxin -inf inf];
elseif isempty(boxin)
    box = [-inf inf -inf inf];
end
if isempty(fit_data)
    fit_leg = '';
end     
if ~isa(smooths_init,'function_handle')
    smooths_init = @(x,y) y;
end
if ~isa(dcs,'function_handle')
    smooths_init = @(x,y) y;
end
nfig = max_fig_num();
plot = plotop1(1);

if (length(plotop1)==2)
    plotsingle = plotop1(2);
    plotsinglesingle = 0;
elseif (length(plotop1)==3)
    plotsingle = plotop1(2);
    plotsinglesingle = plotop1(3);
else
    plotsingle = 0;
    plotsinglesingle = 0;
end
if ~iscell(plotop2)
    plotop2 = num2cell(plotop2);
end
[plotreg, plotrefsub, plotdcsub, plotinteg] = plotop2{:}; % [plot as is, plot ref substracted, plot dc substracted, plot integral]

ns = size(scells,2);
[subrow, subcol] = subplot_min_rectangle(ns);
allleg = cell(1,ns);
allx = NaN(1, ns);
ally = NaN(1, ns);
alldiffdcx = NaN(1, ns);
alldiffdcy = NaN(1, ns);
allintegdiffdcy = NaN(1, ns);

ssize = 0;
if ref
    if ~isempty(trimleg)
        refleg = scells{1,end}(trimleg);
    else
        refleg = scells{1,end};
    end
    reflegsm = [refleg ' smoothed'];
    sref = scells{2,end};
    dvec = ~excludedata(sref(:,1),sref(:,2),'box',box);
    refx = sref(dvec,1);
    refy = sref(dvec,2);
    if any(imag(refy)>0)
        refy = horzcat(real(refy), imag(refy));
    end
    refysm = smooths_init(refx, refy);
    lref = size(refysm, 1);
    dcy = dcs(refx, refysm);
    alldiffleg = cell(1,ns-ref);
    alldiffdcleg = cell(1,ns-ref);
    alldiffx = NaN(1, ns-ref);
    alldiffy = NaN(1, ns-ref);
    
end
hpf = {};
hpp = {};
for i = 1:ns
    if ~isempty(trimleg)
        sleg = scells{1,i}(trimleg);
    else
        sleg = scells{1,i};
    end
    s = scells{2,i};
    if ~isempty(headerparam) && headerparam(1)>0
        hpf = scells{3,i}(1:headerparam,2);
    end
    if ~isempty(headerparam) && length(headerparam)==2 && headerparam(2)>0
        hpp = scells{3,i}(1:headerparam,2);
    end
    dvec = ~excludedata(s(:,1),s(:,2),'box',box);
    sx = s(dvec,1);
    sy = s(dvec,2);
    if any(imag(sy)>0)
        sy = horzcat(real(sy), imag(sy));
    end
    if ~isempty(fit_data)
        [fit_leg, fit_func, fit_plot_range] = fit_data(sx, sy, hpf{:});
    end
    ls = size(sx, 1);
    ssize = max(ssize, ls);
    if ref&&(i<ns)
        if ~isempty(fit_data)
            ref_fit_leg = fit_data(refx, refy, hpf{:});
        end   
        if plotreg
            plotx = horzcat(sx, refx);
            ploty = horzcat(sy, refy);
            plotleg = {[sleg,newline,fit_leg], [refleg,newline,ref_fit_leg]};
            if plotsingle
                ax = line_plot(plotx, ploty, hpp{:}, plotleg, [nfig+1 subrow subcol i]);
                if ~isempty(fit_data)
                    hold(ax(1),'on');
                    fplot(fit_func, fit_plot_range,'Parent',ax(1));
                    hold(ax(1),'off');
                end
            end
            if plotsinglesingle
                ax = line_plot(plotx, ploty, hpp{:}, plotleg, nfig + 1001);
                if ~isempty(fit_data)
                    hold(ax(1),'on');
                    fplot(fit_func, fit_plot_range,'Parent',ax(1));
                    hold(ax(1),'off');
                end
            end
        end
        if plotrefsub||plotdcsub
            l = min(length(refx), ls);
            diffx = NaN(l, 1);
            diffy = NaN(l, 1);
            diffdcy = NaN(l, 1);
            meanxdiff = abs(mean(sx(2:end)-sx(1:end-1)));
            k = 0;
            for j=1:ls
                xdiff = abs(refx - sx(j));
                if min(xdiff)<meanxdiff
                    k = k+1;
                    diffx(k) = sx(j);
                    diffy(k) = sy(j) - refysm(j);
                    diffdcy(k) = sy(j) - dcy(j);
                end
            end
            diffleg = [sleg,' refsm substracted',newline,fit_data(diffx, diffy, hpf{:})];
            if plotrefsub
                if plotsingle
                    line_plot(diffx, diffy, hpp{:}, diffleg, [nfig+2 subrow subcol i])
                end
                if plotsinglesingle
                    line_plot(diffx, diffy, hpp{:}, diffleg, nfig + 1002)
                end
            end
            alldiffleg{i} = diffleg;
            alldiffx(end+1:ssize,:) = NaN;
            alldiffy(end+1:ssize,:) = NaN;
            ldiff = length(diffx);
            alldiffx(1:ldiff,i) = diffx;
            alldiffy(1:ldiff,i) = diffy;
        end
     
    else
        if plotrefsub||plotdcsub
            dcy = dcs(sx, sy);
            diffx = sx;
            diffdcy = sy - dcy;
        end
        if plotreg
            regleg = [sleg newline fit_leg];
            if plotsingle
                ax = line_plot(sx, sy, hpp{:}, regleg, [nfig+1 subrow subcol i]);
                if ~isempty(fit_data)
                    hold(ax(1),'on');
                    fplot(fit_func, fit_plot_range,'Parent',ax(1));
                    hold(ax(1),'off');
                end
            end
            if plotsinglesingle
                ax = line_plot(sx, sy, hpp{:}, regleg, nfig + 1000 + i);
                if ~isempty(fit_data)
                    hold(ax(1),'on');
                    fplot(fit_func, fit_plot_range,'Parent',ax(1));
                    hold(ax(1),'off');
                end
            end
        end
    end
    
    if plotrefsub||plotdcsub
        diffdcleg = [sleg,' DC substracted',newline,fit_data(diffx, diffdcy, hpf{:})];
        alldiffdcleg{i} = diffdcleg;
        alldiffdcx(end+1:ssize,:) = NaN;
        alldiffdcy(end+1:ssize,:) = NaN;
        ldiff = length(diffx);
        alldiffdcx(1:ldiff,i) = diffx;
        alldiffdcy(1:ldiff,i) = diffdcy;
    end
    
    allintegdiffdcy(end+1:ssize,:) = NaN;
    if plotinteg
        if plotdcsub
            integdiffdcy = vec_integ(diffx, diffdcy);
        else
            integdiffdcy = vec_integ(sx, sy);
        end
        allintegdiffdcy(1:ldiff-1,i) = integdiffdcy;
    end
    
    if plotdcsub 
        if plotsingle
            line_plot(diffx, diffdcy, hpp{:}, diffdcleg, [nfig+3 subrow subcol i])
        end
        if plotsinglesingle
            line_plot(diffx, diffdcy, hpp{:}, diffdcleg, nfig + 1003)
        end
    end
    if plotinteg 
        if plotsingle
            integplot(diffx(1:end-1), integdiffdcy, sleg, [nfig+4 subrow subcol i])
        end
        if plotsinglesingle
            integplot(diffx(1:end-1), integdiffdcy, sleg, nfig + 1004)
        end
    end
    
    allleg{i} = [sleg newline fit_leg];
    allx(end+1:ssize,:) = NaN;
    ally(end+1:ssize,:) = NaN;
    allx(1:ls,i) = sx;
    ally(1:ls,i) = sy(:,1);
        
end

if ref
    allx(1:lref,end+1:end+2) = refx*[1 1];
    ally(1:lref,end+1) = refysm;
    allleg{end+1} = reflegsm;
    ally(1:lref,end+1) = dcy;
    allleg{end+1} = 'DC Signal';
end
if plot
    line_plot(allx, ally, allleg, nfig+5);
    if ref&&plotrefsub
    	line_plot(alldiffx, alldiffy, alldiffleg, nfig+6)
    end
    if plotdcsub
        line_plot(alldiffdcx, alldiffdcy, allleg, nfig+7)
    end
    if plotinteg
        integplot(alldiffdcx, allintegdiffdcy, allleg, nfig+8)
    end
end
s_diffdc = {allleg(1:size(alldiffdcx, 2)) alldiffdcx alldiffdcy};
s_diffdc_integ = {allleg(1:size(alldiffdcx, 2)) alldiffdcx allintegdiffdcy};
end

