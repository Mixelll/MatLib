% fileprop = {'txt', '	', [0 0]}; %{'format', 'delimeter', [row_shift column_shift OPTIONAL_to_row OPTIONAL_to_column]}

CV_path = ['C:\Users\' getenv('USERNAME') '\Google Drive\EE MSc\200413 gr-Si Schottky b1\2T Left from center\all CV\try'];
Ct_path = 'C:\Users\admin2\Google Drive\MATLAB Drive\vars\B5 b7 150um 4\B5 b7 150um 4 DLTS1Up';
Ct_path = 'D:\Google Drive\MATLAB Drive\vars\B5 b7 150um 4\B5 b7 150um 4 80-300 DLTS up V0-1';
% Ct_path = 'D:\Google Drive\MATLAB Drive\vars\B5 b7 150um 4\B5 b7 150um 4 DLTS1Up';
Ct_path = 'D:\Google Drive\MATLAB Drive\vars\B5 b7 150um 4\B5 b7 150um 4 80-300 DLTS up PT1';
Ct_path = 'C:\Users\admin2\Google Drive\MATLAB Drive\vars\IB2 CL 100um 10\IB2 CL 100um 10 80-300 DLTS up';

DAQ = true; 
Plotter = false;

if DAQ
    TimeZone = 2; % Hours
    StartTime = 0; % Seconds
    EndTime = Inf; % Seconds
    DataSelect.demods.sample_auxin0_avg = true; DataSelect.imps.sample_param0_avg = true; DataSelect.imps.sample_param1_avg = true;
    [CtCellsFrames, tDelta, T_vs_t] = MFIA_DAQ_MatchData(Ct_path, TimeZone, StartTime, EndTime, DataSelect);
    s_plot(T_vs_t(:,1)/60,T_vs_t(:,2), '.', {}, 'Sample Temperature vs Time', '', '', 'Time [min]', 'Temperature [K]','','','',1, 0)
else
Index = ''; % index (int or str): for example if your files end with a number *ind* and some string *str* after it. input in formet 'indstr' e.g. '5)'. write 'm' instead of ind
SortBy = ''; % sortby (array): sort the files by the numbers that are located (it searches) in the positions you enter in sortby
Save =  0;
paraminheader = {'T'};
CtCells0 = data_folder_read(Ct_path, SortBy, Index, paraminheader, Save); %(path, {type, delimeter, RC}, sortby, ind)
end

%% Frame Segmenter
if DAQ
    [DataTableS,DataCellsStruct,TransitionTable, TransTableSummarized, SetTable] = MFIA_DAQ_FrameSegmenter(CtCellsFrames);
end

%% Global Fit, Plot and Save options
Device = 'B5 b7 150um 4';
% Device = 'IB2 CL 100um 10';
Title = [Device ' C-t'];
SetName = ''; % only applies to un-structured data
FontS.label = 16;
FontS.legend = 16;
PlotOpt = [0 0 0 10]; % [plot, plot single in one window, plot single in multiple windows, plot every #/Param]
PlotHidden = false;
PlotType = @plot;
ScreenSize = [100 100 1000 1000];
LPH = [];
trimheader = [];
FitProp = {'Robust','Bisquare'};
PlotPropCt = {'.'};
PlotPropRes = {'.', 'MarkerSize',20};
FitPlotProp = {'vLineWidth',2};
SaveNameFuncCellLegend = @(f) {findobj(f, 'Type','Line').DisplayName};
SaveNameFuncCellTitle = @(f) {findobj(f, 'Type','Axes', 'UserData',1).Title.String}; 
SaveFormat = {'jpeg', 'fig'};
lblT = 'Temperature [K]';
lblt = 'Time [sec]';
lblC = 'Capacitance [F]';
lblE = 'Emission Time [sec]';
lblN = 'Concentration [cm^{-3}]';
lblS = 'DLTS Signal [a.u.]';
LegNE = 'legend(''Location'',''northeast'')';
LegNW = 'legend(''Location'',''northwest'')';
C0Frac = 100;
% Physical constants
kB = 8.617e-5;

%% Frames Plot
if DAQ
    SignalSlect = {'Capacitance'};
    Range = [65*tDelta inf];
    Range = [];
    TemperatureRange = [0 300.6 1];
    TemperatureRange = [];
%     TemperatureRange = [220 230 1];
    SetNumbers = [20:10:175];
    SetNumbers = [];
    [DataTableFP, tcTemp] = MFIA_DAQ_Frames(CtCellsFrames, 'DataSelect',SignalSlect, 'Range',Range, 'TemperatureRange',TemperatureRange, 'SetNumbers',SetNumbers);
    TitleAddTemp = '';
end
PlotRangeTemp = [65*tcTemp inf];
PlotRangeTemp = [0 inf];
SparsePlot = 15;
rplot = @(sx, sy, leg, subp) s_plot(sx, sy, PlotPropCt, {}, [Title TitleAddTemp], leg, subp, lblt, lblC, FontS, PlotType, PlotHidden, ScreenSize, 0); % s_plot(sx1, sy1, plotproperties, charttitle, legendd, subplott, xlbl, ylbl, fontsize, varargin)
Figures = plot_charts2(DataTableFP, rplot, [1 SparsePlot], '', '', '' ,PlotRangeTemp, 'SpareParameter','T'); % (ScellsIn, PlotFunc, PlotOpt, FitFunc, PlotFit, FitPlotProp, BoxIn, varargin)

%% Segments Plot
if DAQ
    Range = [100*tDelta inf];
    CurveSelect = {["F0T-1"]};
%     CurveSelect = {"F-0.3T-1"};
%     CurveSelect = {["F-0.5T-1"]};
%     CurveSelect = {@(Fm0p3Tm1,Fm0p5Tm1) Fm0p3Tm1-Fm0p5Tm1};
%     CurveSelect = {@(F0Tm1,Fm0p3Tm1) F0Tm1-Fm0p3Tm1};
    TemperatureRange = [200 203 1];
%     TemperatureRange = [220 230 1];
    SetNumbers = [20:10:175];
    SetNumbers = [];
    MovMeanTemp = 1;
    [DataTableSegmentsTemp, CtCellsTemp, tcTemp] = MFIA_DAQ_Process(DataTableS, 'Range',Range, 'CurveSelect',CurveSelect, 'TemperatureRange',TemperatureRange, 'MovMean', MovMeanTemp, 'SetNumbers',SetNumbers);
    TitleAddTemp = ' -';
    for c = CurveSelect
        cc = c{:};
        if isstring(cc)
            TitleAddTemp = [TitleAddTemp ' ' convertStringsToChars(cc)];
        elseif isa(cc, 'function_handle')
            fnstr = func2str(cc);
            TitleAddTemp = [TitleAddTemp ' ' fnstr(find(fnstr==')',1)+1:end)];
        elseif ischar(cc)
            TitleAddTemp = [TitleAddTemp ' ' cc];
        end
    end
end
PlotRangeTemp = [65*tcTemp inf];
PlotRangeTemp = [0 inf];
rplot = @(sx, sy, leg, subp) s_plot(sx, sy, PlotPropCt, {}, [Title TitleAddTemp ', MovAvg=' num2str(MovMeanTemp)], leg, subp, lblt, lblC, FontS, PlotType, PlotHidden, ScreenSize); % s_plot(sx1, sy1, plotproperties, charttitle, legendd, subplott, xlbl, ylbl, fontsize, varargin)
[Figures] = plot_charts2(CtCellsTemp, rplot, 1, '', '', '' ,PlotRangeTemp); % plot_charts2(scells, line_plot, plotop1, trimleg, fit_data, headerparam, boxin)

%% Data processing
MovMean = 1;
if Plotter
N_signals = 3;
SignalSlect = {'Capacitance'};
SegmentSelect = {'low'};
DataRng = [];
    [CtCells, CtCellsStruct, tc]= MFIA_Plotter_Process(CtCells0,N_signals,SignalSlect,SegmentSelect,DataRng,MovMean);
end
if DAQ
    SetNumbers = [11:1:1000];
    Range = [30*tDelta inf];
%     Range = [0 inf];
    TemperatureRange = [0 300.6 1];
%     TemperatureRange = [250 inf];
%     CurveSelect = {["F0T-1"]};
%     CurveSelect = {"F-0.3T-1"};
    CurveSelect = {["F-0.5T-1"]};
%     CurveSelect = {@(F0Tm1,Fm0p05Tm1) F0Tm1-Fm0p05Tm1};
%     CurveSelect = {@(Fm0p3Tm1,Fm0p5Tm1) Fm0p3Tm1-Fm0p5Tm1};
%     CurveSelect = {@(F0Tm1,Fm0p3Tm1) F0Tm1-Fm0p3Tm1};
    [DataTableSP, CtCells, tc] = MFIA_DAQ_Process(DataTableS, 'CurveSelect',CurveSelect, 'TemperatureRange',TemperatureRange, 'Range',Range, 'MovMean',MovMean, 'SetNumbers',SetNumbers);
end


%% Extract DLTS signal
%General options
DLTS_Analyze = 0; % 1 - run analyze and signal plot, 1 - run only signal plot
DLTS_SAVE = 1; % save graphs
% Signal extract and plot options
CtPlotOnce = 1; % plot C-t curves once
PlotFit = 0; % Plot exp(-t/tau) fit for each C-t curve

S_PlotEach = 0; % plot DLTS Signal curves for each rate window USUALLY FALSE
 
N_Signals_Plot = 8; % Select # of rate windows for which a DLTS signal will be plotted 

SignalMode = 0; % dC=C(t1)-C(t2). 0 - dC, 1 - dC/C(end), 2 - dC/(C(0)-C(end)), 3 - 2*N * dC/(C(end)*(exp(t1/tau)-exp(t2/tau)))
DecayUp = 1; % use [] or '' to return to default
SurfaceStates = 1; % Calculate surface trap concentration instead of volume concentration
InterfaceWidth = 0;
S_MovMean = 1;
% Analyze options
CcVec = [0.1 0.15 0.2 0.25 0.3 0.4 0.5]; % T_peak calc by top Cc% of Signal height. low Cc->only high points. you'll get one graph for each Cc_i in Cc vec    x     for each 01 selection below
T_weight_plot = 1; % plot selected Signal-T points that serve (abs(Signal)) as weight points in TpeAK = mean(T, weights=abs(Signal>Cc))
T_peak_plot = 1; % plot temperatures of signal peaks.
e1_fit_plot = 1; % Plot and fit to extract activation energy (E) and pre-exp (A) factor - lin fit [log(e1/T^2) vs 1/T]; e1= A*T^2 * exp(-E/kBT)

Time = now; DTSTR = datestr(Time, 'yyyy-mm-dd HH-MM');
YLim_S = {'ylim([-0.1 1.6])'}; YLim_S = {};


N=1.3e17;
RateWindowCell = {};

switch SignalMode
    case 0
        SignalModeStr = ' Absolute: C(t1)-C(t2)';
        lblS = 'DLTS Signal [F]';
    case 1
        SignalModeStr = ' Relative: C(t1)-C(t2)/C(\infty)';
    case 2
        SignalModeStr = ' Normed: C(t1)-C(t2)/\DeltaC(0)';
    case 3
        SignalModeStr = ' Trap cm^{-3}: 2*N * C(t1)-C(t2)/C(\infty)/e^{t_1}-e^{t_2}';
end

CtFNames = fieldnames(CtCells);
NCtFNames = numel(CtFNames);

if CtPlotOnce
    dplot = @(sx, sy, leg, subp) s_plot(sx, sy, PlotPropCt, Commands, [Title ', Ct\_MovAvg=' num2str(MovMean)], leg, subp, lblt, lblC, FontS, PlotType, PlotHidden, ScreenSize);% s_plot(sx1, sy1, plotproperties, charttitle, legendd, subplott, xlbl, ylbl, fontsize, varargin)
    Figures = plot_charts2(CtCells, dplot, [1 0 0 15], [], [], [] ,[], 'TrimLeg',trimheader, 'HeaderParam',LPH, 'SpareParameter','T'); % plot_charts2(scells, line_plot, plotop1, trimleg, fit_data, headerparam, boxin)
    if DLTS_SAVE
        if isstruct(CtCells)
            ii = floor(numel(Figures)/NCtFNames);
            for i = 1:length(CtFNames)
                SaveDirOuter = [Ct_path '\' [regexprep(CtFNames{i}, {'m','p'}, {'-','\.'}) ' '] ', SigMode=' num2str(SignalMode) ', Ct_MovAvg=' num2str(MovMean) ', S_MovAvg=' num2str(S_MovMean) ' ' DTSTR];
                SaveFigures(Figures((i-1)*ii+1:min(i*ii,numel(Figures))), SaveDirOuter, {}, SaveFormat, ['C-t MovAvg=' num2str(MovMean)]); 
            end
        else
            SaveDirOuter = [Ct_path '\' SetName ', SigMode=' num2str(SignalMode) ', Ct_MovAvg=' num2str(MovMean) ', S_MovAvg=' num2str(S_MovMean) ' ' DTSTR];
            SaveFigures(Figures, SaveDirOuter, {}, SaveFormat, ['C-t MovAvg=' num2str(MovMean)]); 
        end
        
    end 
end
% RateWindowCell{end+1} = [100*tc 0.5]';
% RateWindowCell{end+1} = [[1*tc 20*tc]' [10*tc 200*tc]' [20*tc 400*tc]' [50*tc 1000*tc]' [100*tc 2000*tc]' [150*tc 3000*tc]' [200*tc 4000*tc]'];
% RateWindowCell{end+1} = [[65*tc 650*tc]' [80*tc 800*tc]' [100*tc 1000*tc]' [150*tc 1500*tc]' [200*tc 2000*tc]' [250*tc 2500*tc]' [300*tc 3000*tc]'  [350*tc 3500*tc]' [450*tc 4500*tc]'];
% RateWindowCell{end+1} = [[3*tc 9*tc]' [10*tc 30*tc]' [30*tc 90*tc]' [50*tc 150*tc]' [100*tc 300*tc]' [200*tc 600*tc]' [500*tc 1500*tc]' [800*tc 2400*tc]' [1100*tc 3300*tc]' [1500*tc 4500*tc]'];
% RateWindowCell{end+1} = [[15*tc 75*tc]' [30*tc 150*tc]' [60*tc 300*tc]' [80*tc 400*tc]' [100*tc 500*tc]' [150*tc 750*tc]' [200*tc 1000*tc]' [250*tc 1250*tc]' [300*tc 1500*tc]' [400*tc 2000*tc]' [600*tc 3000*tc]' [800*tc 4000*tc]'];
RateWindowCell{end+1} = RateWindowArray(5, 30*tc, 9400*tc, 100, @LogSpace);
RateWindowCell{end+1} = RateWindowArray(10, 30*tc, 9400*tc, 100, @LogSpace);
RateWindowCell{end+1} = RateWindowArray(20, 30*tc, 9400*tc, 100, @LogSpace);

PlotRange = [0 0.36];
for RWC = RateWindowCell
RateWindow = RWC{:};
Tall = {};
Sall = {};
PeakLoc = {};
% PeakLoc = zeros(1,size(RateWindow,2));
Legall = cell(size(RateWindow,2),NCtFNames);
FileAddall = cell(size(RateWindow,2),NCtFNames);
i = 0;
Ratio = mean(RateWindow(2,:)./RateWindow(1,:));
for RW_i = RateWindow
i = i+1;
% PlotRange = RW_i;
xSp = PlotRange(2)-PlotRange(1);
Commands = {['xlim([' num2str(roundn(PlotRange(1) -0.05*xSp, floor(log10(xSp)-1))) ' ' num2str(roundn(PlotRange(2) +0.2*xSp, floor(log10(xSp)-1))) ']);'], 'legend(''Location'',''southeast'')'};

DLTS = @(x,y) DLTS_Signal(x, y, RW_i, N, mean(y(end-round(length(y)/C0Frac):end)), SignalMode, 'MinLegend',true, 'DecayUp',DecayUp, 'SurfaceStates',SurfaceStates, 'InterfaceWidth',InterfaceWidth);

dplot = @(sx, sy, leg, subp) s_plot(sx, sy, PlotPropCt, Commands, [Title ' - RateWindow= ' sprintf('%0.2g ',RW_i) ', DLTS Signal Mode - ' SignalModeStr ', Ct_MovAvg=' num2str(MovMean)], leg, subp, lblt, lblC, FontS, PlotType, PlotHidden, ScreenSize);% s_plot(sx1, sy1, plotproperties, charttitle, legendd, subplott, xlbl, ylbl, fontsize, varargin)

[Figures, Ct_DLTS_Out] = plot_charts2(CtCells, dplot, PlotOpt, DLTS, PlotFit, FitPlotProp ,PlotRange, 'TrimLeg',trimheader, 'HeaderParam',LPH, 'SpareParameter','T'); % plot_charts2(scells, line_plot, plotop1, trimleg, fit_data, headerparam, boxin)
j = 0;
for cc = fn_struct2cell(Ct_DLTS_Out)
    j = j+1;
    DLTS_Cells = cc{2};
    Ct_FieldName = cc{4};

T = [];
S = [];
Tmin = 0;
Tmax = 310;
for c = DLTS_Cells
    Ti = c{3}{2};
    if Tmin<= Ti && Ti<=Tmax 
        T(end+1) = Ti;
        coeff = c{4};
        S(end+1) = coeff{2,strcmpi(coeff(1,:),'S')};
    end
end
S = movmean(S,S_MovMean);
% Sc=S;
% Sc(abs(S) <= abs(Cc)) = 0;

FontS = FontS;
FontS.legend = 16;

TSp = T(end)-T(1);
cmd = ['xlim([' num2str(roundn(T(1) -0.05*TSp, floor(log10(TSp)-1))) ' ' num2str(roundn(T(end) +0.05*TSp, floor(log10(TSp)-1))) ']);'];
cmd1p = {LegNE cmd};

if ~isempty(Ct_FieldName)
    TitleAdd = [' - ' regexprep(Ct_FieldName, {'m','p'}, {'-','\.'})];
    FileAdd = [regexprep(Ct_FieldName, {'m','p'}, {'-','\.'}) ' '];
else
    TitleAdd = SetName;
    FileAdd = SetName;
end
FileAddall{i,j} = FileAdd;
Leg = ['RateWindow = ' sprintf('%0.2g ',RW_i)];
Legall{i,j} = Leg;
if S_PlotEach
    [~, Figures(end+1)] = s_plot(T, S, PlotPropRes, [{LegNE} cmd YLim_S], [Device TitleAdd ' - DLTS Signal' SignalModeStr ', Ct\_MovAvg=' num2str(MovMean) ', S\_MovAvg=' num2str(S_MovMean)], Leg, [], lblT, lblS, FontS, PlotType, false, true);
end
if j>length(Tall) || isempty(Tall{j})
    Tall{j} = T.';
    Sall{j} = S.';
%     PeakLoc{j} = sum(T.*Sc)/sum(Sc);
else
    Tall{j} = [Tall{j} T.'];
    Sall{j} = [Sall{j} S.'];
%     PeakLoc{j} = [PeakLoc{j} sum(T.*Sc)/sum(Sc)];
end

SaveDir = [Ct_path '\' FileAdd ', SigMode=' num2str(SignalMode) ', Ct_MovAvg=' num2str(MovMean) ', S_MovAvg=' num2str(S_MovMean) ' ' DTSTR '\Ratio=' num2str(Ratio,2)  '\RateWindow= ' sprintf('%0.2g ',RW_i)];
if DLTS_SAVE && numel(Figures), SaveFigures(Figures, SaveDir, SaveNameFuncCellLegend, SaveFormat); end
% close all
end
end

% FigureS = gobjects(0);
% FigureDLTS = gobjects(0);

for k = 1:j
    SaveDirOuter = [Ct_path '\' FileAddall{1,k}  ', SigMode=' num2str(SignalMode) ', Ct_MovAvg=' num2str(MovMean) ', S_MovAvg=' num2str(S_MovMean) ' ' DTSTR '\Ratio=' num2str(Ratio,2) ];
    if N_Signals_Plot
        ii = ceil(size(Tall{k},2)/N_Signals_Plot);
        [AX, FigureS] = s_plot(Tall{k}(:,1:ii:end), Sall{k}(:,1:ii:end), PlotPropRes, [{cmd} YLim_S], [Device TitleAdd ' - DLTS Signal' SignalModeStr ', Ratio=' num2str(Ratio,2) ', Ct\_MovAvg=' num2str(MovMean) ', S\_MovAvg=' num2str(S_MovMean)], Legall(1:ii:end,k), [], lblT, lblS, FontS, PlotType, false, true);
    end
    if DLTS_SAVE, SaveFigures(FigureS, SaveDirOuter, {}, SaveFormat, ['DLTS Signal Mode ' num2str(SignalMode)]); end
end
%% Analyze DLTS signal
if DLTS_Analyze

% Cc = 0.75;
T_RangeVec_peak = [[0 ;270]];
T_RangeVec_fit = [[0 ;400] [0 ;210] [210 ;400]];
T_RangeVec_fit = [[0 ;400]];

for Cc = CcVec
for T_Range_peak = T_RangeVec_peak
for T_Range_fit = T_RangeVec_fit
T_Range_str = sprintf('%d ',T_Range_fit);
SaveDirAllInner = [SaveDirOuter '\T Weighted by DLTS Peak=' num2str(Cc*100) '%, T_Range=' T_Range_str(1:end-1)];
for k = 1:j
    T_mat =  Tall{k};
    T = mean(T_mat,2); TL = T>=T_Range_peak(1) & T<=T_Range_peak(2);
    T_mat = T_mat(TL,:);
    RW = RateWindow;
    S_mat = Sall{k}(TL,:);
    S_mat_max = max(abs(S_mat));
    Sc=S_mat;
    LMat = abs(S_mat) <= abs((1-Cc)*S_mat_max);
    Sc(LMat) = NaN;
    if T_weight_plot
    ii = ceil(size(Tall{k},2)/N_Signals_Plot);
    if N_Signals_Plot
        [~, FigureDLTS] = s_plot(T_mat(:,1:ii:end), Sc(:,1:ii:end), PlotPropRes, [{cmd} YLim_S], [Device TitleAdd ' - DLTS Signal' SignalModeStr ', Ratio=' num2str(Ratio,2) ', Ct\_MovAvg=' num2str(MovMean) ', S\_MovAvg=' num2str(S_MovMean) ' - Peak=' num2str(Cc*100) '% of Signal'], Legall(1:ii:end,k), [], lblT, lblS, FontS, PlotType, false, true);
        if DLTS_SAVE, SaveFigures(FigureDLTS, SaveDirAllInner, {}, SaveFormat, 'T Peak Weight Points'); end 
    end
    end
    Sc(LMat) = 0;
    Tpeak = sum(T_mat.*Sc)./sum(Sc);
    dRW = RW(1,:)-RW(2,:);
    if T_peak_plot, [~, FigureDLTS] = s_plot(dRW, Tpeak, 'ro', {}, [Device TitleAdd ' - DLTS Signal' SignalModeStr ', Ratio=' num2str(Ratio,2) ', Ct\_MovAvg=' num2str(MovMean) ', S\_MovAvg=' num2str(S_MovMean) ' - Peak T vs \Deltat  - Peak=' num2str(Cc*100) '%'], '', [], '\Deltat [sec]', lblT, FontS, PlotType, false, true);
    if DLTS_SAVE, SaveFigures(FigureDLTS, SaveDirAllInner, {}, SaveFormat, 'Peak T'); end 
    end
    if e1_fit_plot
    e1 = log(RW(1,:)/RW(2,:))./dRW;
    T1 = 1000./Tpeak;
    [lbl, fun, span, Fcell] = DLTS_Fit(Tpeak,e1, 'Range',T_Range_fit); E = Fcell{2,1}; A = Fcell{2,2};
    [AX, FigureDLTS, ttl, ~, ~, ylbl] = s_plot(T1,  e1./Tpeak.^2, 'ro', LegNE, [Device TitleAdd ' - DLTS Signal' SignalModeStr ', Ratio=' num2str(Ratio,2) ', Ct\_MovAvg=' num2str(MovMean) ', S\_MovAvg=' num2str(S_MovMean) ' - e_{1}/T^{2} vs 1000/T - Peak=' num2str(Cc*100) '%'], ['Measured Data' newline lbl], [], '1000/T [K^{-1}]', '$\frac{e_{1}}{T^{2}}\mathrm{\hspace{1mm}[sec^{-1}]}$', FontS, @semilogy, false, true);
    ylbl.Interpreter = 'latex'; ylbl.FontSize = 25; %ttl.Interpreter = 'latex'; 
    hold(AX(1), 'on')
    fplot(fun, span)
    if DLTS_SAVE, SaveFigures(FigureDLTS, SaveDirAllInner, {}, SaveFormat, 'log e1 by T^2 vs 1000 by T'); save([SaveDirOuter '\Data and Fit'],'Tpeak','e1','E','A'); end
    end
end
end
end
end
% hold off
end
DLTS_Analyze;
end

%% C-t Fit (tau and N)
Time = now;
YLim_tau1 = {'ylim([-0.018 0])'};
YLim_Nt1 = {'ylim([-0.018 0])'};
YLim_tau2 = {'ylim([-0.018 0])'};
YLim_Nt2 = {'ylim([-0.018 0])'};
YLim_tau1 = {}; YLim_Nt1 = {}; YLim_tau2 = {}; YLim_Nt2 = {};
PlotFit = true;
N=1.3e17;
NT = [1];

FitRange = [[10*tc 2.1]' [100*tc 2.1]' [1000*tc 2.1]' [10000*tc 2.1]' [10*tc 100*tc]' [10*tc 1000*tc]' [10*tc 10000*tc]' [100*tc 1000*tc]'];
FitRange = [[0 1000*tc]' [300*tc 3000*tc]' [1000*tc 10000*tc]' [2400*tc 24000*tc]'];
FitRange = [[60*tc 200*tc]' [200*tc 600*tc]' [500*tc 1500*tc]' [1000*tc 3000*tc]' [2000*tc 6000*tc]' [4000*tc 12000*tc]' [8000*tc 24000*tc]'];
FitRange = [[60*tc 1400*tc]' [200*tc 4000*tc]' [500*tc 10000*tc]' [800*tc 16000*tc]' [1200*tc 24000*tc]'];
FitRange = [100*tc 0.46]';
PlotRange = [0 0.46];

for NT_i = NT
Tall = {};
Nt1all = {};
tau1all = {};
Nt2all = {};
tau2all = {};
Legall = cell(0,size(FitRange,2));
FileAddall = cell(0,size(FitRange,2));
i = 0;
for FT_i = FitRange
i = i+1;
PlotRange = FT_i;
xSp = PlotRange(2)-PlotRange(1);
Commands = {['xlim([' num2str(roundn(PlotRange(1) -0.05*xSp, floor(log10(xSp)-1))) ' ' num2str(roundn(PlotRange(2) +0.3*xSp, floor(log10(xSp)-1))) ']);'], 'legend(''Location'',''southeast'')'};

Ct = @(x,y) Ct_Nt_tau_fit(x,y,FT_i,mean(y(end-round(length(y)/C0Frac):end)),N,NT_i,FitProp);
 
rplot = @(sx, sy, leg, subp) s_plot(sx, sy, PlotPropCt, Commands, Title, leg, subp, lblt, lblC, FontS, PlotType, PlotHidden, ScreenSize);% s_plot(sx1, sy1, plotproperties, charttitle, legendd, subplott, xlbl, ylbl, fontsize, varargin)

[Figures, CtFitOut] = plot_charts2(CtCells, rplot, PlotOpt, Ct, PlotFit, FitPlotProp,PlotRange, 'TrimLeg',trimheader, 'HeaderParam',LPH, 'SpareParameter','T'); % plot_charts2(scells, line_plot, plotop1, trimleg, fit_data, headerparam, boxin)
j = 0;
for cc = fn_struct2cell(CtFitOut)
    j = j+1;
    Ct_FittedCells = cc{2};
    Ct_FieldName = cc{4};

T = [];
Nt1 = [];
tau1 = [];
if NT_i>1
    Nt2 = [];
    tau2 = [];
end
Tmin = 0;
Tmax = 310;
for c = Ct_FittedCells
    Ti = c{3}{2};
    if Tmin<= Ti && Ti<=Tmax 
        T(end+1) = Ti;
        coeff = c{4};
        Nt1(end+1) = coeff{2,strcmpi(coeff(1,:),'Nt1')};
        tau1(end+1) = coeff{2,strcmpi(coeff(1,:),'tau1')};
        if NT_i>1
            Nt2(end+1) = coeff{2,strcmpi(coeff(1,:),'Nt2')};
            tau2(end+1) = coeff{2,strcmpi(coeff(1,:),'tau2')};
        end
    end
end
TSp = T(end)-T(1);
cmd = ['xlim([' num2str(roundn(T(1) -0.05*TSp, floor(log10(TSp)-1))) ' ' num2str(roundn(T(end) +0.05*TSp, floor(log10(TSp)-1))) ']);'];
cmd1p = {LegNE cmd};
cmd2p = {LegNW cmd ; LegNE cmd};
% [~, Figures(end+1)] = s_plot(T, tau1, {'.', 'MarkerSize',40}, 'legend(''Location'',''northeast'')', [Title ' - Trap Emission Times and Trap Concentration'], {'\tau_{1}'}, [], 'Temperature [K]', 'Emission Time [sec]', Font, PlotType, false, true);
if ~isempty(Ct_FieldName)
    TitleAdd = [' - ' regexprep(Ct_FieldName, {'m','p'}, {'-','\.'})];
    FileAdd = [regexprep(Ct_FieldName, {'m','p'}, {'-','\.'}) ', '];
else
    TitleAdd = SetName;
    FileAdd = SetName;
end
FileAddall{i,j} = FileAdd;
Leg = ['FitRange = ' sprintf('%0.2g ',FT_i)];
Legall{i,j} = Leg;
[~, Figures(end+1)] = s_plot(T, tau1, PlotPropRes, cmd2p, [Device TitleAdd ' FitRange= ' sprintf('%0.2g ',FT_i) '- Trap 1 Emission Times and Concentration'], {'\tau_{1}',0,'Nt_{1}'}, [], lblT, {lblE lblN}, FontS, PlotType, false, true, '',T, Nt1);
if NT_i==2
    [~, Figures(end+1)] = s_plot(T, tau2, PlotPropRes, cmd2p, [Device TitleAdd ' FitRange= ' sprintf('%0.2g ',FT_i) '- Trap 2 Emission Times and Concentration'], {'\tau_{2}',0,'Nt_{2}'}, [], lblT, {lblE lblN}, FontS, PlotType, false, true, '',T, Nt2);
end

if j>length(Tall) || isempty(Tall{j})
    Tall{j} = T.';
    tau1all{j} = tau1.';
    Nt1all{j} = Nt1.';
    if NT_i==2
        tau2all{j} = tau2.';
        Nt2all{j} = Nt2.';
    end
else
    Tall{j} = [Tall{j} T.'];
    tau1all{j} = [tau1all{j} tau1.'];
    Nt1all{j} = [Nt1all{j} Nt1.'];
    if NT_i==2
        tau2all{j} = [tau1all{j} tau2.'];
        Nt2all{j} = [Nt1all{j} Nt2.'];
    end
end

SaveDir = [Ct_path '\' FileAdd 'C-t Fit, MovAvg=' num2str(MovMean) ' ' datestr(Time, 'yyyy-mm-dd HH-MM') '\NT='  num2str(NT_i) '\FitRange= ' sprintf('%0.2g ',FT_i) ', MovAvg=' num2str(MovMean)]; 
SaveFigures(Figures, SaveDir, SaveNameFuncCellLegend, SaveFormat)
end
end
for k = 1:j
    L = Legall(:,k);
    TwoSetLegend = @(a,b) [cellfun(@(c) strcat(a,c), L, 'UniformOutput', false)' {0} cellfun(@(c) strcat(b,c), L, 'UniformOutput', false)'];
    SaveDirOuter = [Ct_path '\' FileAddall{1,k}  'C-t Fit, MovAvg=' num2str(MovMean) ' ' datestr(Time, 'yyyy-mm-dd HH-MM') '\NT='  num2str(NT_i)];
    [~, FiguresCtFit] = s_plot(Tall{k}, tau1all{k}, PlotPropRes, [cmd1p YLim_tau1], [Device TitleAdd ' - C-t Fit exp' num2str(NT_i) ' - tau1'], L, [], lblT, lblE, FontS, PlotType, false, true);
    [~, FiguresCtFit(end+1)] = s_plot(Tall{k}, Nt1all{k}, PlotPropRes, [cmd1p YLim_Nt1], [Device TitleAdd ' - C-t Fit exp' num2str(NT_i) ' - Nt1'], L, [], lblT, lblN, FontS, PlotType, false, true);
    [~, FiguresCtFit(end+1)] = s_plot(Tall{k}, tau1all{k}, PlotPropRes, [cmd2p [YLim_tau1; YLim_Nt1]], [Device TitleAdd ' - C-t Fit exp' num2str(NT_i) ' - tau1, Nt1'], TwoSetLegend('\tau_{1}, ','Nt_{1}, '), [], lblT, {lblE lblN}, FontS, PlotType, false, true, '', T, Nt1all{k});

    if NT_i==2
        [~, FiguresCtFit(end+1)] = s_plot(Tall{k}, tau2all{k}, PlotPropRes, [cmd1p YLim_tau2], [Device TitleAdd ' - C-t Fit exp2 - tau2'], L, [], lblT, lblE, FontS, PlotType, false, true);
        [~, FiguresCtFit(end+1)] = s_plot(Tall{k}, Nt2all{k}, PlotPropRes, [cmd1p YLim_Nt2], [Device TitleAdd ' - C-t Fit exp2 - Nt2'], L, [], lblT, lblN, FontS, PlotType, false, true);
        [~, FiguresCtFit(end+1)] = s_plot(Tall{k}, tau2all{k}, PlotPropRes, [cmd2p [YLim_tau2; YLim_Nt2]], [Device TitleAdd ' - C-t Fit exp' num2str(NT_i) ' - tau2, Nt2'], TwoSetLegend('\tau_{2}, ','Nt_{2}, '), [], lblT, {lblE lblN}, FontS, PlotType, false, true, '', T, Nt2all{k});
        [~, FiguresCtFit(end+1)] = s_plot(Tall{k}, tau1all{k}, PlotPropRes, [cmd2p [YLim_tau1; YLim_tau2]], [Device TitleAdd ' - C-t Fit exp' num2str(NT_i) ' - tau1, tau2'], TwoSetLegend('\tau_{1}, ','\tau_{2}, '), [], lblT, {lblE lblE}, FontS, PlotType, false, true, '', T, tau2all{k});
        [~, FiguresCtFit(end+1)] = s_plot(Tall{k}, Nt1all{k}, PlotPropRes, [cmd2p [YLim_Nt1; YLim_Nt2]], [Device TitleAdd ' - C-t Fit exp' num2str(NT_i) ' - Nt1, Nt2'], TwoSetLegend('Nt_{1}, ','Nt_{2}, '), [], lblT, {lblN lblN}, FontS, PlotType, false, true, '', T, Nt2all{k});
    end

    SaveFigures(FiguresCtFit, SaveDirOuter, SaveNameFuncCellTitle, SaveFormat, 'Regexp',' - ', 'After',true)
end

end

