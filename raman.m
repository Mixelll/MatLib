function [] = raman()
%path, index, title, ref, r, plotop,
fileprop = {'txt', '', [42 0]}; %{'format', 'delimeter', [row_shift column_shift OPTIONAL_to_row OPTIONAL_to_column]}
DCns = [10 3]; %dcnsegments nmins
user_def_DC = 1000; % place a string here to use a signal-generated DC
user = 'z1x38';
user = 'sz1x3';
%user = 'Mixel';
%path = ['C:\Users\' user '\Google Drive\EE MSc\In with the New\Measurements\RAMAN\Ilya Group\Michael\Implantation 1 SOI\3'];
path = ['C:\Users\' user '\Google Drive\EE MSc\In with the New\Manufacturing\Implantations\Si\Raman\First Implantation SOI\3'];
%path = ['C:\Users\' user '\Google Drive\EE MSc\In with the New\Manufacturing\Implantations\Si\Raman\Second Implantation Si'];
index = 3; % '' if not index, 'm' for max 
title = 'Raman Spectra of 100keV Ar+ Implanted Si';
%title = '';
plotop1 = [1 0 1]; % [plot, plot single in one window, plot single in multiple windows]
plotop2 = [0 0 0 0]; % [plot as is, plot ref substracted, plot dc substracted, plot integral]
smoothmethod = @smooth;
acqxP = 15*1;
r = 50;
r = r*acqxP;
range = [100 200 400 500 500 600 100 700];
[Rcells, ref] = data_folder_read(path, fileprop{:}, index);
rsmooth = @(sx, sy) smooth_s(sx, sy, r, smoothmethod, DCns(1));
dcs = @(sx, sy) dc_sig(sx, sy, r, DCns, user_def_DC); % do not input user_def_DC to use a signal-generated DC


rplot = @(sx, sy, leg, subp) s_plot(sx, sy, title, leg, subp,'Shift [1/cm]', 'Intensity [Counts]', 15); % (sx, sy, charttitle, leg, subp, xlbl, ylbl, fontsize)

rintegplot = @(sx, sy, leg, subp) s_plot(sx, sy, ['Integrated' title], leg, subp, 'Shift [1/cm]', 'Integrated Intensity [Counts/cm]', 15);
[DiffDC] = plot_charts(Rcells, rplot, rintegplot, rsmooth, dcs, ref, plotop1, plotop2, ''); % (scells, s_plot, integplot, smooths_init, dcs, ref, plotop1, plotop2, fit)

[Signalarea, SignalareaScaled, SignalareaScaledFractioned] = signal_area(DiffDC{end-2}, DiffDC{end-1}, DiffDC{end}, range, acqxP, '[Counts]*[1/cm]')
NoiseScaled = s_noise(DiffDC{end-2}, DiffDC{end-1}, DiffDC{end}, rsmooth, range, acqxP, '[Counts]*[1/cm]')

end

