fileprop = {'txt', '	', [0 0]}; %{'format', 'delimeter', [row_shift column_shift OPTIONAL_to_row OPTIONAL_to_column]}
id = '4';
path = ['E:\Gr-Si sch Implanted\310825 IB4 5E14 150C N2 annealed\center right\pairs before and after 100C\' id];
%path = vladpath;
index = '';
sortby = '';
saveop =  0;
title = ['IV Gr-Si scht 1E14 imp before and after 100C N2 anneal - center right ' id];
font = 25;
plotop1 = [1 1]; % [plot, plot single in one window, plot single in multiple windows]
plotop2 = [1 0 0 0]; % [plot as is, plot ref substracted, plot dc substracted, plot integral]


paraminheader = {};
lph = length(paraminheader);
trimheader = '';
[VIcells, ref] = data_folder_read(path, fileprop{:}, sortby, index, paraminheader, saveop); %(path, {type, delimeter, RC}, sortby, ind)
% sortby (array): sort the files by the numbers that are located (it searches) in the positions you enter in sortby
% ind (int or str): for example if your files end with a number *ind* and some string *str* after it. input in formet 'indstr' e.g. '5)'. write 'm' instead of ind

% INPUT ALL VALUES  in CM "
syms N vfb Ar VB n V a n0
n0=1e12;
%a='';
N=1;
n0=1;
n=1;
f=1;    
R = 100e-4;
A = pi*R^2;
q = 1.6e-19;
e0 = 8.854e-14;
es = 11.68;

T = 296;




%          dV/dln(I) = Rs=I + n*kT
plot_range = [0 1]; %V
segments1 = [0 0.15 0.3 0.45 0.6 0.75 0.9 1.05 1.2 1.35 1.5]; % [fit_range1_i fit_range1_f fit_range2_i fit_range2_f...]
segments1 = [0 1];
segments2 = [segments1(1:end-1);segments1(2:end)];
fit_range_global = [0.5 1];
r_n_vs_segments = false;
title1 = ['Rs extraction'];
title2 = ['Ideality factor extraction'];


Fit_outcells_Rs_n = cell(size(VIcells)+[1 0]);
Fit_outcells_Rs_n(1,1:end) = VIcells(1,1:end);
for r = segments2
    fit_range = r';
    fit_range = fit_range_global;
%     plot_range = fit_range;
    fit_model_R = @(x,y) I_Rs_n_schot(x,y,fit_range,T,1);
    fit_model_n = @(x,y) I_Rs_n_schot(x,y,fit_range,T,2);
    fit_model_1 = @(x,y) I_Rs_n_schot(x,y,fit_range,T); 
    fit_model_2 = @(x,y) I_I0_n_schot_2(x,y,fit_range,fit_model_R(x,y),T);
    rplot_1 = @(sx, sy, leg, subp) s_plot(sy(1:end-1) + diff(sy)/2, diff(sx)./diff(log(sy)),'', title1, {leg, 0, 'V=data(I)','V=data(I)-Rs*I'}, subp, 'Current [I]', {'dV/dln(I)','Voltage [V]'}, font,'',sy,sx,sy,sx-sy*fit_model_R(sx,sy)); % (sx, sy, charttitle, leg, subp, xlbl, ylbl, fontsize)
    plot_charts(VIcells, rplot_1, '', '', '', ref, plotop1, plotop2, trimheader, fit_model_1, [lph lph], plot_range); % (scells, s_plot, integplot, smooths_init, dcs, ref, plotop1, plotop2, fit, plot_range)
    rplot_2 = @(sx, sy, leg, subp) s_plot(sx-sy*fit_model_R(sx,sy), log(sy./(1-exp(-(sx-sy*fit_model_R(sx,sy))/0.0257))), '', title2, leg, subp, 'Voltage -IRs [V]', 'ln(I/(1-exp(-(V-IR)/kT))', font); % (sx, sy, charttitle, leg, subp, xlbl, ylbl, fontsize)
    plot_charts(VIcells, rplot_2, '', '', '', ref, plotop1, plotop2, trimheader, fit_model_2, [lph lph], plot_range); % (scells, s_plot, integplot, smooths_init, dcs, ref, plotop1, plotop2, fit)
    for c = 1:size(Fit_outcells_Rs_n,2)
        onedeviceVI = VIcells{2,c};
        V = onedeviceVI(1:end,1);
        I = onedeviceVI(1:end,2);
        Fit_outcells_Rs_n{2,c} = vertcat(Fit_outcells_Rs_n{2,c},[mean([r(1) r(2)]) fit_model_R(V,I)]);
        Fit_outcells_Rs_n{3,c} = vertcat(Fit_outcells_Rs_n{3,c},[mean([r(1) r(2)]) fit_model_n(V,I)]);
    end
end
if r_n_vs_segments
    plotproperties = {'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'b'};
    rplot = @(sx, sy, leg, subp) s_plot(sx, sy, plotproperties, 'Rs vs V', leg, subp, 'Voltage [V]', 'R [ohm]', font);
    plot_charts(Fit_outcells_Rs_n([1 2],:), rplot, '', '', '', ref, plotop1, plotop2, '', '', '', '');
    rplot = @(sx, sy, leg, subp) s_plot(sx, sy, plotproperties, 'n vs V', leg, subp, 'Voltage [V]', '\eta', font);
    plot_charts(Fit_outcells_Rs_n([1 3],:), rplot, '', '', '', ref, plotop1, plotop2, '', '', '', '');
end

%          ln(I/(1-exp(-V/kT)) = ln(I0) + kT/n*V
%fit_range = [0.1 0.25];
%plot_range = [0.1 0.25];
% fit_model = @(x,y) I_I0_n_schot_1(x,y,fit_range); 
% rplot = @(sx, sy, leg, subp) s_plot(sx, log(sy./(1-exp(-sx/0.0257))), title, leg, subp, 'Voltage [V]', 'ln(I/(1-exp(-V/kT))', font); % (sx, sy, charttitle, leg, subp, xlbl, ylbl, fontsize)
% plot_charts(VIcells, rplot, '', '', '', ref, plotop1, plotop2, fit_model, plot_range); % (scells, s_plot, integplot, smooths_init, dcs, ref, plotop1, plotop2, fit)

%          ln(I/(1-exp(-(V-IRs)/kT)) = ln(I0) + kT/n*(V-IRs)
%fit_range = [0.1 0.25];
%plot_range = [0.1 0.25];
% fit_model_R = @(x,y) I_Rs_n_schot(x,y,fit_range,1);
% fit_model = @(x,y) I_I0_n_schot_2(x,y,fit_model_R(x,y),fit_range); 
% rplot = @(sx, sy, leg, subp) s_plot(sx-sy*fit_model_R(sx,sy), log(sy./(1-exp(-(sx-sy*fit_model_R(sx,sy))/0.0257))), title, leg, subp, 'Voltage -IRs [V]', 'ln(I/(1-exp(-(V-IR)/kT))', font); % (sx, sy, charttitle, leg, subp, xlbl, ylbl, fontsize)
% plot_charts(VIcells, rplot, '', '', '', ref, plotop1, plotop2, fit_model, plot_range); % (scells, s_plot, integplot, smooths_init, dcs, ref, plotop1, plotop2, fit)
