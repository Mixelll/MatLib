function [Figures, OutSCell] = plot_charts2(ScellsIn, PlotFunc, PlotOpt, FitFunc, PlotFit, FitPlotProp, BoxIn, varargin)
ColorOrder =... 
         [0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840];
Ncolors = size(ColorOrder,1);
Figures = gobjects(0);

CellStrChar = @(s) iscell(s) || isstring(s) || ischar(s);
ChrStr = @(s) ischar(s) || isstring(s);
CellStrCharFn = @(s) iscell(s) || isstring(s) || ischar(s) || isa(s, 'function_handle');
p = inputParser;
p.KeepUnmatched=true;
p.addParameter('TrimLeg', []); % Select indices of legend chars to display
p.addParameter('HeaderParam', {}); % Get parameters from file header or table
p.addParameter('SpareParameter', ''); % Get parameters from file header or table

p.parse(varargin{:});

TL = p.Results.TrimLeg;
HP = p.Results.HeaderParam;
SP = p.Results.SpareParameter;


if size(BoxIn,2)==2
    Box = [BoxIn -inf inf];
elseif size(BoxIn,1)==2
    Box = [BoxIn.' -inf inf];
elseif isempty(BoxIn)
    Box = [-inf inf -inf inf];
end
if isempty(PlotFit)
    PlotFit = true;
end
if isempty(FitFunc)
    FitLeg = '';
end
if isempty(FitPlotProp)
    FitPlotProp = {};
end

Plot = PlotOpt(1);

if length(PlotOpt)==2 && any([0 1]==PlotOpt(2))
    PlotSingle = PlotOpt(2);
    PlotSingleSingle = 0;
elseif length(PlotOpt)==3 && any([0 1]==PlotOpt(3))
    PlotSingle = PlotOpt(2);
    PlotSingleSingle = PlotOpt(3);
else
    PlotSingle = 0;
    PlotSingleSingle = 0;
end
if length(PlotOpt)==4 || (length(PlotOpt)==2 && all([0 1]~=PlotOpt(2))) || (length(PlotOpt)==3 && all([0 1]~=PlotOpt(3)))
    Sparse = PlotOpt(end);
    SparseFlag = true;
else
    Sparse = 1;
    SparseFlag = false;
end

Param = {};
HPF_flag = false;
HPP_flag = false;
SP_flag = false;
if isstring(HP)
    HP = cellstr(HP);
elseif ~iscell(HP)
    HP = {HP};
end
if ~isempty(HP) && ~isempty(HP{1,1})
    HPF_flag = true;
    Param = HP(1,:);
end
if ~isempty(HP) && size(HP,1)==2 && ~isempty(HP{2,1})
    HPP_flag = true;
    Param = [HP(1,:) HP(2,:)];
end
if ~isempty(SP) && SparseFlag
    SP_flag = true;
    if ~any(strcmp(Param, SP))
        Param{end+1} = SP;
    end
end
if istable(ScellsIn)
    ScellsIn = Table2Cell_PC(ScellsIn, Param);
end
SF = false; % StructFlag
if isstruct(ScellsIn)
    FieldNames = fieldnames(ScellsIn).';
    SF = true;
else
    Scells = ScellsIn;
end

f = 1;
while true
    if SF
        Scells = ScellsIn.(FieldNames{f});
    end
    Nfig = max_fig_num();
    NS = size(Scells,2);
    OutBottomCell = cell(1,NS);
    FittedFuncsCell = cell(1,NS);
    FitPlotRange = cell(1,NS);
    AllLeg = cell(1,NS);
    Allx = NaN(1, NS);
    Ally = NaN(1, NS);
    if PlotSingle
        [subrow, subcol] = subplot_min_rectangle(NS);
    end
    
    Ssize = 0;
    HPF = {}; % params that go into fit function
    HPP = {}; % params that go into plot function
    SPP_last = nan;
    for i = 1:NS
        SparsePlot = true;
        if SP_flag
            SparsePlot = false;
            SPP = Scells{3,i}{strcmp(Scells{3,i}(:,1), SP),2};
            if SPP-SPP_last>Sparse || isnan(SPP_last)
                SPP_last = SPP;
                SparsePlot = true;
            end
        elseif SparseFlag
            SparsePlot = ~mod(i,Sparse);
        end
        if ~isempty(TL)
            Sleg = Scells{1,i}(TL);
        else
            Sleg = Scells{1,i};
        end
        S = Scells{2,i};
        if HPF_flag
%             HPF = [];
%             for c = HeaderParam(1,:)
%                 HPF = [HPF Scells{3,i}(strcmp(Scells{3,i}(:,1), c),2)];
%             end
            HPF = cellfun(@(c) Scells{3,i}(strcmp(Scells{3,i}(:,1), c),2), HP(1,:));
        end
        if HPP_flag
            HPP = cellfun(@(c) Scells{3,i}(strcmp(Scells{3,i}(:,1), c),2), HP(2,:));
        end
       
        dvec = ~excludedata(S(:,1),S(:,2),'box',Box);
        Sx = S(dvec,1);
        Sy = S(dvec,2);
        if any(imag(Sy)>0)
            Sy = horzcat(real(Sy), imag(Sy));
        end
        if ~isempty(FitFunc)
            [FitLeg, FittedFuncsCell{i}, FitPlotRange{i}, OutBottomCell{i}] = FitFunc(Sx, Sy, HPF{:});
        end
        ls = size(Sx, 1);
        Ssize = max(Ssize, ls);
        RegLeg = [Sleg newline FitLeg];
        if PlotSingle && SparsePlot
            [ax, figsbp] = PlotFunc(Sx, Sy, HPP{:}, RegLeg, [Nfig+1 subrow subcol i]);
            if ~isempty(FitFunc) && ~isempty(FittedFuncsCell{i}) && PlotFit
                hold(ax(1),'on');
                fplot(FittedFuncsCell{i}, FitPlotRange{i}, FitPlotProp{:}, 'Parent',ax(1));
                hold(ax(1),'off');
            end
        end
        if PlotSingleSingle && SparsePlot
            [ax, fig] = PlotFunc(Sx, Sy, HPP{:}, RegLeg, Nfig + 1000 + i);
            if ~isempty(FitFunc) && ~isempty(FittedFuncsCell{i}) && PlotFit
                hold(ax(1),'on');
                fplot(FittedFuncsCell{i}, FitPlotRange{i}, FitPlotProp{:}, 'Parent',ax(1));
                hold(ax(1),'off');
            end
            Figures(end+1) = fig;
        end
        if  SparsePlot
            AllLeg{i} = [Sleg newline FitLeg];
            Allx(end+1:Ssize,:) = NaN;
            Ally(end+1:Ssize,:) = NaN;
            Allx(1:ls,i) = Sx;
            Ally(1:ls,i) = Sy(:,1);
        end
    end
    if PlotSingle
        Figures = [figsbp Figures];
    end
    AllBoolVec = cellfun(@isempty, AllLeg);
    AllLeg(AllBoolVec) = [];
    Allx(:,AllBoolVec) = [];
    Ally(:,AllBoolVec) = [];
    if ~isempty(FitFunc)
        FittedFuncsCell(AllBoolVec) = [];
    end
    if Plot && (isempty(FitFunc) || ~PlotFit)
        [~, fig] = PlotFunc(Allx, Ally, AllLeg, Nfig+5);
    elseif Plot && PlotFit
        for i = 1:size(Ally,2)
            [ax, fig] = PlotFunc(Allx(:,i), Ally(:,i), HPP{:}, AllLeg{i}, Nfig+5);
            ax(1).Children(1).Color = ColorOrder(mod(i,Ncolors)+1,:);
            hold(ax(1),'on');
            FitPlotR = [max(FitPlotRange{i}(1),BoxIn(1)) min(FitPlotRange{i}(2),BoxIn(2))];
            fplot(FittedFuncsCell{i}, FitPlotR, FitPlotProp{:}, 'Parent',ax(1), 'Color',ColorOrder(mod(i,Ncolors)+1,:));
        end
        if size(Ally,2)
            hold(ax(1),'off');
        end
    end
    if Plot && any(~AllBoolVec)
        Figures = [fig Figures];
    end
    OutSCell = [Scells ; OutBottomCell];
    if SF
        OutStructSCell.(FieldNames{f}) = OutSCell;
        f = f+1;
        if f>length(FieldNames)
            break
        end
    else
        break
    end
end
if SF
    OutSCell = OutStructSCell;
end
end

