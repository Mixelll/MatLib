function [] = IV()

%path, index, title, ref, r, plotop,

fileprop = {'txt', '	', [0 0]}; %{'format', 'delimeter', [row_shift column_shift OPTIONAL_to_row OPTIONAL_to_column]}
id = '9';
path = ['E:\Gr-Si sch Implanted\310825 IB4 5E14 150C N2 annealed\center right\pairs before and after 100C\' id];
%path = vladpath;
index = '';
sortby = '';
saveop =  0;
title = ['IV Gr-Si scht 1E14 imp before and after 100C N2 anneal - center right ' id];
font = 20;
plotop1 = [1 0]; % [plot, plot single in one window, plot single in multiple windows]
plotop2 = [1 0 0 0]; % [plot as is, plot ref substracted, plot dc substracted, plot integral]

paraminheader = {};
lph = length(paraminheader);
trimglegend = '';

[VIcells, ref] = data_folder_read(path, fileprop{:}, sortby, index, paraminheader, saveop); %(path, {type, delimeter, RC}, sortby, ind, tparam, saveop)
% sortby (array): sort the files by the numbers that are located (it searches) in the positions you enter in sortby
% ind (int or str): for example if your files end with a number *ind* and some string *str* after it. input in formet 'indstr' e.g. '5)'. write 'm' instead of ind
% saveop (int): 1 to save txt collection in one cell array, 2 to save as a matrix for each txt

fit = @(x,y) str_poly_fit(x,y,1,1e-1); %(x,y,poly degree,neglect coefficients that contribute less fraction)
fit = '';
dcs = @(sx, sy) dc_sig(sx, sy, r, DCns);
dcs = '';
range = [-inf inf];
% 
%
plotproperties = {};
rplot = @(sx, sy, leg, subp) s_plot(sx, log10(abs(sy)), plotproperties, title, leg, subp, 'Voltage (V)', 'Current (log(A))', font); % (sx, sy, charttitle, leg, subp, xlbl, ylbl, fontsize)
plot_charts(VIcells, rplot, '', '', dcs, ref, plotop1, plotop2, trimglegend, fit, lph, range);

end




