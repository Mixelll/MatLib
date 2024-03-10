function [] = IV_gr_si_schot_reverse()


% fileprop = {'mat', '	', [0 0]}; %{'format', 'delimeter', [row_shift column_shift OPTIONAL_to_row OPTIONAL_to_column]}
FileType = 'mat';
path = 'C:\Users\user\My Drive\MATLAB Drive\vars\B5 b5 150um 14\B5 b5 150um 14 Vacuum1\DC\IV.mat';
% vladpath = 'C:\Users\admin2\Google Drive\EE MSc\In with the New\ppt\Vlad';
%path = vladpath;
index = '';
sortby = '';
saveop =  0;
title = 'Graphene-Si Schottky IV';
font = 13;
plotop1 = [0 1]; % [plot, plot single in one window, plot single in multiple windows]
plotop2 = [1 0 0 0]; % [plot as is, plot ref substracted, plot dc substracted, plot integral]


if FileType=='mat'
    VIcells = load(path);
    VIcells_fields = fieldnames(VIcells);
    VIcells = VIcells.(VIcells_fields{1});
    VIcells = [VIcells.sweep.V', VIcells.sweep.I'];
    ref = 0;
else
    [VIcells, ref] = data_folder_read(path, sortby, index, 0, saveop, 'FileType', FileType); %(path, sortby, ind, tparam, saveop)
% sortby (array): sort the files by the numbers that are located (it searches) in the positions you enter in sortby
% ind (int or str): for example if your files end with a number *ind* and some string *str* after it. input in formet 'indstr' e.g. '5)'. write 'm' instead of ind
end
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
Ac30 = 132*200;
Ac50 = 150*240;
Ac70 = 168*280;
Ac100 = 196*340;
Ac = Ac100*1e-8 -A/2;
q = 1.6e-19;
e0 = 8.854e-14;
es = 11.68;

range = [-3 0];
fit_model_Rev = @(x,y) I_reverse_a(x,-y,range,f,A,Ar,es,N,VB,a,n0);
%rplot = @(sx, sy, leg, subp) s_plot(sx, sy,'', title, {leg, 0, 'I=data(V)','I=data(V-Rs*I)'}, subp, 'Voltage [V]', 'Current [I]', font,sy,sx-sy*fit_model_Rev(sx,sy));
rplot = @(sx, sy, leg, subp) s_plot(sx, -sy,'', title, leg, subp, 'Voltage [V]', 'Current [I]', font);
plot_charts(VIcells, rplot, '', '', '', ref, plotop1, plotop2, fit_model_Rev, range);

end




