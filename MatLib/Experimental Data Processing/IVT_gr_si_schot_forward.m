function [] = IVT_gr_si_schot_forward()


fileprop = {'txt', '	', [0 0]}; %{'format', 'delimeter', [row_shift column_shift OPTIONAL_to_row OPTIONAL_to_column]}
path = 'C:\Users\admin2\Google Drive\EE MSc\In with the New\Meas\Probing\T C-t\4 - Gr-Si B1 6 Left Center\Vacuum\T Down\IV\All T no hyst';
index = '';
sortby = '';
saveop =  0;
title = 'Graphene-Si Schottky IV';
font = 13;
plotop = [0 1 0]; % [plot, plot single in one window, plot single in multiple windows]


paraminheader = {'T'};
lph = length(paraminheader);
trimheader = [1:6];
[VITcells] = data_folder_read(path, fileprop{:}, sortby, index, paraminheader, saveop); %(path, {type, delimeter, RC}, sortby, ind, tparam, saveop)
% sortby (array): sort the files by the numbers that are located (it searches) in the positions you enter in sortby
% ind (int or str): for example if your files end with a number *ind* and some string *str* after it. input in formet 'indstr' e.g. '5)'. write 'm' instead of ind
% saveop (int): 1 to save txt collection in one cell array, 2 to save as a matrix for each txt

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
k = 8.617e-5;


%          dV/dln(I) = Rs=I + n*kT
plot_range = [0 1]; %V
fit_range_R_n = [0.8 1]; %V
fit_range_I0_n = [0.05 1]; %V
fit_range_richard = [120 301]; %T


% plot_range = fit_range_R_n;
fit_model_R = @(x,y,T) I_Rs_n_schot(x,y,fit_range_R_n,T,1);
fit_model_n = @(x,y,T) I_Rs_n_schot(x,y,fit_range_R_n,T,2);
fit_model_I0 = @(x,y,T) I_I0_n_schot_2(x,y,fit_range_I0_n,fit_model_R(x,y,T),T,1);
plot_fit_model_R = @(x,y,T) I_Rs_n_schot(x,y,fit_range_R_n,T);
plot_fit_model_I0_R = @(x,y,T) I_I0_n_schot_2(x,y,fit_range_I0_n,fit_model_R(x,y,T),T);
plot_fit_model_richard = @(x,y) I0_T_Richardson(x,y,fit_range_richard,A);
% We get Rs and n from here
rplot_1 = @(sx, sy, T, leg, subp) s_plot(sy(1:end-1) + diff(sy)/2, diff(sx)./diff(log(sy)),'', title, {leg, 0, 'V=data(I)','V=data(I)-Rs*I'}, subp, 'Current [I]', {'dV/dln(I)','Voltage [V]'}, font,'',sy,sx,sy,sx-sy*fit_model_R(sx,sy,T)); % (sx, sy, charttitle, leg, subp, xlbl, ylbl, fontsize)
% plot_charts2(VITcells, rplot_1, plotop2, trimheader, plot_fit_model_R, [lph lph], plot_range); % (scells, line_plot, plotop1, trimleg, fit_data, headerparam, boxin)
% We get I0 and n from here
rplot_2 = @(sx, sy, T, leg, subp) s_plot(sx-sy*fit_model_R(sx,sy,T), log(sy./(1-exp(-(sx-sy*fit_model_R(sx,sy, T))/(k*T)))), '', title, leg, subp, 'Vapplied-IRs [V]', 'ln(I/(1-exp(-(V-IR)/kT))', font); % (sx, sy, charttitle, leg, subp, xlbl, ylbl, fontsize)
% plot_charts2(VITcells, rplot_2, plotop2, trimheader, plot_fit_model_I0_R, [lph lph], plot_range); % (scells, line_plot, plotop1, trimleg, fit_data, headerparam, boxin)

Fit_outcells_Rs_n_lnI0byT2 = cell([2 3]);
Fit_outcells_Rs_n_lnI0byT2(1,1:end) = {'6 Left center'};
for c = 1:size(VITcells,2)
    oneTVI = VITcells{2,c};
    T = VITcells{3,c}{1,2};
    V = oneTVI(1:end,1);
    I = oneTVI(1:end,2);
    Fit_outcells_Rs_n_lnI0byT2{2,1} = vertcat(Fit_outcells_Rs_n_lnI0byT2{2,1},[T fit_model_R(V,I,T)]);
    Fit_outcells_Rs_n_lnI0byT2{2,2} = vertcat(Fit_outcells_Rs_n_lnI0byT2{2,2},[T fit_model_n(V,I,T)]);
    Fit_outcells_Rs_n_lnI0byT2{2,3} = vertcat(Fit_outcells_Rs_n_lnI0byT2{2,3},[T fit_model_I0(V,I,T)]);
end

plotproperties = {'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'b'};
rplot = @(sx, sy, leg, subp) s_plot(sx, sy, plotproperties, 'Rs vs T', leg, subp, 'Temperature [K]', 'R [ohm]', font);
% plot_charts2(Fit_outcells_Rs_n_lnI0byT2(:,1), rplot, plotop, '', '', '', '');
rplot = @(sx, sy, leg, subp) s_plot(sx, sy, plotproperties, 'n vs T', leg, subp, 'Temperature [K]', 'Ideality Factor \eta', font);
% plot_charts2(Fit_outcells_Rs_n_lnI0byT2(:,2), rplot, plotop, '', '', '', '');
rplot = @(sx, sy, leg, subp) s_plot(1./sx, log(sy./sx.^2), plotproperties, 'ln(I0/T^2) vs 1/T', leg, subp, 'Temperature^{-1} [K^{-1}]', 'ln(I0/T^2)', font);
plot_charts2(Fit_outcells_Rs_n_lnI0byT2(:,3), rplot, plotop, '', plot_fit_model_richard, '', '');

% Trushin
% limits = {'VB',[0.1 1.5 2],'n',[1 5 100]};
% limits = [limits {'a',[-2 0.1 2]} ];
% limits = [limits {'n0',[1e10 1e13 1e14]}];
% plot_fit_model_Trushin = @(x,y) I_gr_si_schot_forward_Trushin(x,y,fit_range_Trushin,f,T,fit_model_R(x,y),A,VB,n,n0,outfittedvar_Trushin,limits{:});
% rplot_3 = @(sx, sy, leg, subp) s_plot(sx-sy*fit_model_R(sx,sy), sy, '', title, leg, subp, 'Vapplied-IRs [V]', 'Current [I]', font); % (sx, sy, charttitle, leg, subp, xlbl, ylbl, fontsize)
% plot_charts(VITcells, rplot_3, '', '', '', ref, plotop1, plotop2, plot_fit_model_Trushin, plot_range); % (scells, s_plot, integplot, smooths_init, dcs, ref, plotop1, plotop2, fit)





end
