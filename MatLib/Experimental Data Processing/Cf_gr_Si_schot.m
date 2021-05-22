fileprop_txt_tabdelimeted = {'txt', '	', [0 0]}; %{'format', 'delimeter', [row_shift column_shift OPTIONAL_to_row OPTIONAL_to_column]}
fileprop_MFIA = {'txt', ';', [5 0]}; %{'format', 'delimeter', [row_shift column_shift OPTIONAL_to_row OPTIONAL_to_column]}
fileprop = fileprop_MFIA;
path = 'C:\Users\admin2\Google Drive\EE MSc\In with the New\Meas\Probing\Gr-Si not implanted\MFIA\Plot';
index = ''; % ind (int or str): for example if your files end with a number *ind* and some string *str* after it. input in formet 'indstr' e.g. '5)'. write 'm' instead of ind
sortby = ''; % sortby (array): sort the files by the numbers that are located (it searches) in the positions you enter in sortby
saveop =  0;    
title = 'Graphene-Si Schottky B5 150um C-f';
font = 15;
plotop1 = [1 0]; % [plot, plot single in one window, plot single in multiple windows]
plotop2 = [1 0 0 0]; % [plot as is, plot ref substracted, plot dc substracted, plot integral]

fit = [];
paraminheader = {};
lph = length(paraminheader);
trimheader = [];
range = [100.5 inf]; % in data units
plot_type = @semilogx; % @insert_desirable_function

[Cfcells, ref] = data_folder_read(path, fileprop{:}, sortby, index, paraminheader, saveop); %(path, {type, delimeter, RC}, sortby, index, paraminheader, saveop)

rplot = @(sx, sy, leg, subp) s_plot(sx, sy, '', title, leg, subp, 'Frequency [Hz]', 'Capacitance [F]', font, plot_type);% s_plot(sx1, sy1, plotproperties, charttitle, legendd, subplott, xlbl, ylbl, fontsize, plot_type, varargin)

plot_charts(Cfcells, rplot, 0, '', '', ref, plotop1, plotop2, trimheader, fit, lph, range); % plot_charts(scells, line_plot, integplot, smooths_init, dcs, ref, plotop1, plotop2, fit_data, headerparam, boxin)





