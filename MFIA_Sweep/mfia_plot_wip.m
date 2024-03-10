% Define files to load (struct 3D format) and check if exists
FolderPath = ['C:\Users\' getenv('USERNAME') '\My Drive\MATLAB Drive\vars'];
% FolderPath = 'D:\Google Drive\MATLAB Drive\vars';
% EXAMPLE NameCells = {'Sample_Name' {'Measurement_ID'}}; EXAMPLE
NameCells = {'B5 b5 150um 4' {'A'} 'B5 b5 150um 9' {'A' 'DLCP' 'DLCP1' 'M 296K'} 'B5 b5 150um 14' {'B' 'C' 'M 296K' 'M 340K' 'Vacuum1'} ... % {SampleName1 {SweepID1 ... SweepIDn} ... SampleNameN {SweepID1 ... SweepIDn}}
    'B5 b7 150um 4' {'300KVacuum1' '300KVacuum2' 'DLCP2T1' 'DLCP4T1' 'M 170K' 'M 200K' 'M 230K' 'M 260K' 'M 298K' 'M2T1' 'M2T2' 'M4T1' 'M4T2'} 'B5 b7 150um 9' {'300KVacuum1' 'A' 'M 170K' 'M 200K' 'M 230K' 'M 260K' 'M 298K' 'M 340K'}...
    'B5 b7 150um 14' {'300KVacuum1'} 'B5 b7 150um 19' {'300KVacuum1'} 'IB2 BL 100um 2' {'A'} 'IB2 BL 100um 3' {'A'} 'IB2 CL 100um 10' {'300KVacuum1' 'DLCP2T1' 'DLCP4T Forward1' 'DLCP4T1' 'DLCP4T2'}...
    'IB2 CR 100um 2' {'A' 'DLCP4T1' 'DLCP4T2' 'M1'} 'IB3 5E14 H-annealed TL 100um 2' {'A' 'B' 'DLCP1'}};
NameCells = {'B5 b7 150um 4' {'M4T1'}};
NameCells
NameCellsIn = reshape(NameCells, [2 numel(NameCells)/2]);
NameCellsIn2 = {};
for cn = NameCellsIn
    for cnn = cn{2}
        NameCellsIn2(:,end+1) = {cn{1}; cnn(1)};
        FileName = [cn{1} ' ' cnn{1}];
        RealPath = [FolderPath '\' cn{1} '\' FileName '\' FileName];
        if ~isfile([RealPath '.mat'])
            error(['File not found at path' newline RealPath '.mat'])
        else
            disp(['File found at path' newline RealPath '.mat'])
        end    
    end
end

% Process (define range and limits) and optionally plot data
% Device Area
A = (150^2 * pi)*1e-8;
% Save data slice plots
SaveDataPlots = false;
% Select data to plot (leave empty to plot all data in 3D struct) and define plot aliases (from IA parametric model)
plt_select_data = {};
% plt_select_data(:,end+1) = {'param0'; 'Resistance [ohm]'};
plt_select_data(:,end+1) = {'param1'; 'Capacitance [F]'};
% plt_select_data(:,end+1) = {'absz'; 'Impedance Magnitude [ohm]'};
% plt_select_data(:,end+1) = {'phasez'; 'Impedance Phase'};
% Select axis range - data outside range is discarded from here onward (applies to all subsequent fits)
ax_range = {};
ax_range(:,end+1) = {'frequency'; [0 inf]};
ax_range(:,end+1) = {'amplitude'; [0 inf]}; 
ax_range(:,end+1) = {'offset'; [-inf inf]};
% Select value range - data outside range is replaced with NaN from here onward (applies to all subsequent fits)
val_range = {};
val_range(:,end+1) = {'param0'; [0 inf]};
val_range(:,end+1) = {'param1'; [0 inf]};
% Define slice planes
slice_planes = {};
slice_planes(:,end+1) = {'frequency'; [1e2 1e3 1e4 1e5 5e5]};
slice_planes(:,end+1) = {'amplitude'; []}; 
slice_planes(:,end+1) = {'offset'; [0.1 0.2 0.5 1 1.5 0 -0.1 -0.2 -0.5 -1 -1.5]};
% Add plot formatting commands
% all: plt_cmds(2,end+1) =  {'comand'}; '' for ' inside string.
% target: plt_cmds(:,end+1) =  {'target' ; 'comand'};
plt_cmds = {};
% plt_cmds(2,end+1) = {'grid on'};
plt_cmds(2,end+1) = {'colorbar(''eastoutside'')'};
plt_cmds(2,end+1) = {'colormap(s,''turbo'')'};
plt_cmds(2,end+1) = {'colormap(s,interp1(colormap(s), 1:(length(colormap(s))/998):length(colormap(s))))'};
plt_cmds(:,end+1) = {'param0'; 'colormap(s,[compress_array_exp(colormap(s),10,10, ''Reverse'', true, ''Interp'', true)])'};
% plt_cmds(2,end+1) = {'colormap(s,[0 0 0; colormap(s)])'};
% plt_cmds(2,end+1) = {'colormap([colormap; 1 1 1])'};
plt_cmds(:,end+1) = {'param0'; 'caxis([0 1e7])'};
plt_cmds(:,end+1) = {'param1'; 'caxis([0 inf])'};

plt_log_freq = true; % true for log plot frequency
% Load File
for cn = NameCellsIn2
SampleName = cn{1};
SweepID = cn{2}{1};
FileName = [SampleName ' ' SweepID];
SetName = '3D_sweep_order';
RealPath = [FolderPath '\' SampleName '\' FileName '\' FileName];
%RealPath = '' % Or put the path yourself

LoadedFile = load(RealPath)
fNames = fieldnames(LoadedFile);
DataStruct0 = LoadedFile.(fNames{contains(fNames, SetName)});
Order = DataStruct0.order;
% StrBeforeOrder = 'order_'; % input the string just before the first parameter
% Order = MFIA_get_order_from_path(fNames{contains(fNames,SetName)}, StrBeforeOrder)
title = RealPath(sum(find(RealPath=='\', 1, 'last'))+1:end)

SavePlotPath = [FolderPath '\' SampleName '\' title];
% Plot slice planes of 3D data
if plt_log_freq
    Order{contains(Order, 'frequency')} = 'log_frequency';
    ax_range(:,contains(ax_range(1,:), 'frequency')) = {'log_frequency'; log10(ax_range{2,contains(ax_range(1,:), 'frequency')})};
    slice_planes(:,contains(slice_planes(1,:), 'frequency')) = {'log_frequency'; log10(slice_planes{2,contains(slice_planes(1,:), 'frequency')})};
end
try
[PlotFig, SubPlots, DataStruct,~,~,ND] = process_plot_struct_data3D(DataStruct0, Order, slice_planes, 'title', title, 'ax_range', ax_range, 'val_range', val_range, 'plt_select_data', plt_select_data, 'plt_cmds', plt_cmds);
PlotGUI = MFIA_PGUI(SubPlots,ND);
if SaveDataPlots
    StartTime = datestr(now, 'yyyy-mm-dd HH-MM-SS');
    if ~isfolder([SavePlotPath '\Slice Plots'])
        mkdir([SavePlotPath '\Slice Plots']);
    end
    f=figure('Position',get(0,'Screensize'));
    copyobj(get(PlotFig, 'Children'),f)
    f.Visible = 'on';
    saveas(f, [SavePlotPath '\Slice Plots\' FileName ' slice_plot ' StartTime '.fig']);
    f.Visible = 'off';
    saveas(f, [SavePlotPath '\Slice Plots\' FileName ' slice_plot ' StartTime '.bmp']);
    close(f)
    if ~isfolder([SavePlotPath '\GUI'])
        mkdir([SavePlotPath '\GUI']);
    end   
    saveas(PlotGUI, [SavePlotPath '\GUI\' FileName ' GUI ' StartTime '.fig']);
end
catch 
end
end

% Fit data along a certain dimension
% CV fit

% Fit, plot and save options
Save_CV_plot = true;

CV_plot_fit = [];
CV_plot_fit.func = @plot;
CV_plot_fit.visible = 0;
CV_plot_fit.screensize = 0;
CV_plot_fit.savepath = SavePlotPath;
CV_plot_fit.properties = {'Marker','o', 'LineStyle','none'};
% CV_plot_fit = [];

CV_fit_select = {'capacitance'};
CV_fit_axis = 'offset';

CV_val_range = {};
CV_val_range(:,end+1) = {''; [0 inf]};
CV_val_range(:,end+1) = {''; [0 inf]};
% Fit function and options
CV_offset_Range = [-1 -0.5];
CV_str = ['Fit range= ' sprintf('%0.3g ',CV_offset_Range)];
CV_FitProp = {'Robust','Bisquare'};
syms V N Vb n
limits_A = {'Vb',[0 0.5 1],'N',[1e16 2e18 5e18],'n',[1 3 10]};
% limits_A = {};
n=1;
CV_str = [strAdd('n',n, 'N',N, 'Vb',Vb, CV_FitProp{:}) ' ' CV_str];
es = 11.68;
CV_fit_func = @(x,y) C_schot_fit_A(x,y,CV_offset_Range,A,es,N,Vb,n,CV_FitProp,limits_A{:});

[CV_fig, CV_sbp, CV_struct_array, CV_FitAppendPath, CV_FitName] = fit_cell_3D(DataStruct.data, DataStruct.axes, CV_fit_func, CV_fit_axis, 'plt_select_data',CV_fit_select, 'title',title, 'plot_fit',CV_plot_fit, 'val_range',CV_val_range, 'str',CV_str);
CV_GUI = MFIA_PGUI(CV_sbp,2);
if Save_CV_plot
    CV_save_path = [SavePlotPath '\' CV_FitAppendPath];
    mkdir(CV_save_path);
    save([CV_save_path '\CV fit out ' CV_FitName], 'CV_struct_array')
    f=figure('Position',get(0,'Screensize'));
    copyobj(get(CV_fig, 'Children'),f)
    f.Visible = 'on';
    saveas(f, [CV_save_path '\CV fit ' CV_FitName '.fig']);
    f.Visible = 'off';
    saveas(f, [CV_save_path '\CV fit ' CV_FitName '.bmp']);
    close(f)
    saveas(CV_GUI, [CV_save_path '\CV fit GUI ' CV_FitName '.fig']);
end





function str = strAdd(varargin)
str = '';
for i=1:2:length(varargin)
    c = varargin{i+1};
    if isnumeric(c)
        str = [varargin{i} '=' num2str(c,2) ' ' str];
    elseif isstring(c) || ischar(c)
        str = [varargin{i} '=' c ' ' str];
    end
end
end