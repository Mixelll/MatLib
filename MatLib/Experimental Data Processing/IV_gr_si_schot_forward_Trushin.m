function [] = IV_gr_si_schot_forward_Trushin()


fileprop = {'txt', '	', [0 0]}; %{'format', 'delimeter', [row_shift column_shift OPTIONAL_to_row OPTIONAL_to_column]}
path = 'C:\Users\admin2\Google Drive\EE MSc\200413 gr-Si Schottky b1\2T Left from center\all IV\try';
vladpath = 'C:\Users\admin2\Google Drive\EE MSc\In with the New\ppt\Vlad';
%path = vladpath;
index = '';
sortby = '';
saveop =  0;
title = 'Graphene-Si Schottky IV';
font = 13;
plotop1 = [0 1]; % [plot, plot single in one window, plot single in multiple windows]
plotop2 = [1 0 0 0]; % [plot as is, plot ref substracted, plot dc substracted, plot integral]



[VIcells, ref] = data_folder_read(path, fileprop{:}, sortby, index, saveop); %(path, {type, delimeter, RC}, sortby, ind)
% sortby (array): sort the files by the numbers that are located (it searches) in the positions you enter in sortby
% ind (int or str): for example if your files end with a number *ind* and some string *str* after it. input in formet 'indstr' e.g. '5)'. write 'm' instead of ind

% INPUT ALL VALUES  in CM "
syms N Vb Ar VB n V a n0
%n0=1e12;
a=0;
R = 100e-4;
A = pi*R^2;
Ac30 = 132*200;
Ac50 = 150*240;
Ac70 = 168*280;
Ac100 = 196*340;
Ac = Ac100*1e-8 -A/2;
q = 1.6e-19;
e0 = 8.854e-14;
es = 11.68;



%          dV/dln(I) = Rs=I + n*kT
plot_range = [0 1.5]; %V
fit_range_R_n = [0.5 1.5]; %V
fit_range_Trushin = [0.5 1.5]; %V
f=1;
T='';
outfittedvar_Trushin = '';
% plot_range = fit_range_R_n;
fit_model_R = @(x,y) I_Rs_n_schot(x,y,fit_range_R_n,1);
fit_model_n = @(x,y) I_Rs_n_schot(x,y,fit_range_R_n,2);
plot_fit_model_R = @(x,y) I_Rs_n_schot(x,y,fit_range_R_n);
plot_fit_model_I0_R = @(x,y) I_I0_n_schot_2(x,y,fit_model_R(x,y),fit_range_R_n);
% rplot_1 = @(sx, sy, leg, subp) s_plot(sy(1:end-1) + diff(sy)/2, diff(sx)./diff(log(sy)),'', title, {leg, 0, 'V=data(I)','V=data(I)-Rs*I'}, subp, 'Current [I]', {'dV/dln(I)','Voltage [V]'}, font,'',sy,sx,sy,sx-sy*fit_model_R(sx,sy)); % (sx, sy, charttitle, leg, subp, xlbl, ylbl, fontsize)
% plot_charts(VIcells, rplot_1, '', '', '', ref, plotop1, plotop2, plot_fit_model_R, plot_range); % (scells, s_plot, integplot, smooths_init, dcs, ref, plotop1, plotop2, fit, plot_range)
% rplot_2 = @(sx, sy, leg, subp) s_plot(sx-sy*fit_model_R(sx,sy), log(sy./(1-exp(-(sx-sy*fit_model_R(sx,sy))/0.0257))), '', title, leg, subp, 'Vapplied-IRs [V]', 'ln(I/(1-exp(-(V-IR)/kT))', font); % (sx, sy, charttitle, leg, subp, xlbl, ylbl, fontsize)
% plot_charts(VIcells, rplot_2, '', '', '', ref, plotop1, plotop2, plot_fit_model_I0_R, plot_range); % (scells, s_plot, integplot, smooths_init, dcs, ref, plotop1, plotop2, fit)

% Trushin
limits = {'VB',[0.1 1.5 2],'n',[1 5 100]};
limits = [limits {'a',[-2 0.1 2]} ];
limits = [limits {'n0',[1e10 1e13 1e14]}];
plot_fit_model_Trushin = @(x,y) I_gr_si_schot_forward_Trushin(x,y,fit_range_Trushin,f,T,fit_model_R(x,y),A,VB,n,n0,outfittedvar_Trushin,limits{:});
rplot_3 = @(sx, sy, leg, subp) s_plot(sx-sy*fit_model_R(sx,sy), sy, '', title, leg, subp, 'Vapplied-IRs [V]', 'Current [I]', font); % (sx, sy, charttitle, leg, subp, xlbl, ylbl, fontsize)
plot_charts(VIcells, rplot_3, '', '', '', ref, plotop1, plotop2, plot_fit_model_Trushin, plot_range); % (scells, s_plot, integplot, smooths_init, dcs, ref, plotop1, plotop2, fit)





end
