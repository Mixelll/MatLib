function [] = Hall()

%path, index, title, ref, r, plotop,

fileprop = {'txt', '	', [0 0]}; %{'format', 'delimeter', [row_shift column_shift OPTIONAL_to_row OPTIONAL_to_column]}
path = 'C:\Users\sz1x3\Google Drive\EE MSc\In with the New\Measurements\Probing\Graphene Hall bars 050919\plot\Hall\A 2um overlap';
path = 'C:\Users\sz1x3\Google Drive\EE MSc\In with the New\Measurements\Probing\Graphene Hall bars 050919\plot\Hall\B 3um overlap\B 3um overlap recontact';
%path = 'C:\Users\sz1x3\Google Drive\EE MSc\In with the New\Measurements\Probing\Graphene Hall bars 050919\plot\Hall\B 3um overlap\B 3um overlap bot';
index = '';
title = '';
plotop1 = [1 1 0]; % [plot, plot single in one window, plot single in multiple windows]
plotop2 = [1 0 0 0]; % [plot as is, plot ref substracted, plot dc substracted, plot integral]


[Hallcells, ref] = data_folder_read(path, fileprop{:}, index);
fit = @(x,y) str_poly_fit(x,y,1,1e-12);
 
rplot = @(sx, sy, leg, subp) s_plot(sx, sy, title, leg, subp,'Current [A]', 'Voltage [V]', 15); % (sx, sy, charttitle, leg, subp, xlbl, ylbl, fontsize)
plot_charts(Hallcells, rplot, 0, 0, 0, ref, plotop1, plotop2, fit); % (scells, s_plot, integplot, smooths_init, dcs, ref, plotop1, plotop2, plottitle)
% Noise = s_noise(DiffDC{end-2}, DiffDC{end-1}, DiffDC{end}, rsmooth, range, acqxP, '[Counts]*[1/cm]');
% [Signalarea, SignalareaScaled] = signal_area(DiffDC{end-2}, DiffDC{end-1}, DiffDC{end}, range, acqxP, '[Counts]*[1/cm]')

end




