fileprop = {'txt', '	', [0 0]}; %{'format', 'delimeter', [row_shift column_shift OPTIONAL_to_row OPTIONAL_to_column]}
root_path = 'My Drive\MATLAB Drive\vars\';
path = 'D:\Cascade 301\301226 gr-Si Schottky b1 C-V\2T\C-V\100um\try';
%path = 'C:\Users\admin2\Google Drive\EE MSc\200413 gr-Si Schottky b1\2T Left from center\all CV\try';
path = strcat(root_path, 'B5 b5 150um 4\B5 b5 150um 4 A');
index = ''; % ind (int or str): for example if your files end with a number *ind* and some string *str* after it. input in formet 'indstr' e.g. '5)'. write 'm' instead of ind
sortby = ''; % sortby (array): sort the files by the numbers that are located (it searches) in the positions you enter in sortby
saveop =  0;    
title = 'Graphene-Si Schottky CV';
font = 15;
plotop1 = [0 1]; % [plot, plot single in one window, plot single in multiple windows]
plotop2 = [1 0 0 0]; % [plot as is, plot ref substracted, plot dc substracted, plot integral]

paraminheader = {};
lph = length(paraminheader);
trimheader = [];
plot_type = @plot;
[CVcells, ref] = data_folder_read(path, fileprop{:}, sortby, index, paraminheader, saveop); %(path, {type, delimeter, RC}, sortby, index, paraminheader, saveop)

% INPUT ALL VALUES  in CM "
syms N Vb n n0 V
%n0=1e12;
a='';
T='';
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
model = A.*sqrt(q.*es.*e0.*N./(2.*n.*(n.*(Vb-0.0257)-V)));
lin_model = (2.*n.*(n.*(Vb-0.0257)-V))./(q.*es.*e0.*N.*A.^2);


over = true;
outfittedvar = '';

limits_A = {'Vb',[0 0.5 1],'N',[1e17 2e18 5e18],'n',[1 9 10]};
limits_B = [limits_A {'a',[0 0.3 2]} ];
limits_B = [limits_A {'n0',[1e10 1e13 1e14]}];
range = [-2 -1];
f = 1e10;
%for n0=[1e13 1e12 1e11 1e10]
%(x,y,boxin,f,A,es,N,Vb,n,a,n0,varargin)
%CV0 = @(x,y) C_schot_fit_A(x,y,range,f,A,11.68,N,Vb,n,limits_A{:});
%CV1 = @(x,y) C_schot_fit_B1(x,y,range,f,A,11.68,N,Vb,n,a,limits_B{:});
%CV2 = @(x,y) C_schot_fit_B2(x,y,range,f,A,11.68,N,Vb,n,a,limits_B{:});
%CV3 = @(x,y) C_schot_fit_B3(x,y,range,f,A,11.68,N,Vb,n,a,limits_B{:});
%CV4 = @(x,y) C_schot_fit_B4(x,y,range,f,A,11.68,N,Vb,n,a,limits_B{:});

CV0 = @(x,y) C_schot_fit_A(x,y,range,f,A,11.68,N,Vb,n,limits_A{:});
CV1 = @(x,y) C_schot_fit2_B1(x,y,range,f,T,A,11.68,N,Vb,n,a,n0,outfittedvar,limits_B{:});
CV2 = @(x,y) C_schot_fit2_B2(x,y,range,f,T,A,11.68,N,Vb,n,a,n0,outfittedvar,limits_B{:});
CV3 = @(x,y) C_schot_fit2_B3(x,y,range,f,T,A,11.68,N,Vb,n,a,n0,outfittedvar,limits_B{:});
CV4 = @(x,y) C_schot_fit2_B4(x,y,range,f,T,A,11.68,N,Vb,n,a,n0,outfittedvar,limits_B{:});

rplot = @(sx, sy, leg, subp) s_plot(sx, sy, '', title, leg, subp, 'Vbias [V]', 'Capacitance [F]', font, plot_type);% s_plot(sx1, sy1, plotproperties, charttitle, legendd, subplott, xlbl, ylbl, fontsize, plot_type, varargin)
%for CV = {CV1,CV2,CV3,CV4}
for CV = {CV1,CV4}
    plot_charts(CVcells, rplot, 0, '', '', ref, plotop1, plotop2, trimheader, CV{:}, lph, range); % plot_charts(scells, line_plot, integplot, smooths_init, dcs, ref, plotop1, plotop2, fit_data, headerparam, boxin)

end
%end


% f = 1e-20;
% fi = @(x,y) C_schot_fit_1overCsq(x,1./y.^2,range,f,over,A,11.68,N,Vb,n,limits_A{:});
% if over
%     rplot = @(sx, sy, leg, subp) s_plot(sx, sy, title, leg, subp, 'Vbias [V]', 'C^{-2} [F^{-2}]', font);% (sx, sy, charttitle, leg, subp, xlbl, ylbl, fontsize, varargin of LineSpec name-value pairs )
% end
%     
% plot_charts(CVcells, rplot, 0, '', '', ref, plotop1, plotop2, fi); % (scells, s_plot, integplot, smooths_init, dcs, ref, plotop1, plotop2, fit)




