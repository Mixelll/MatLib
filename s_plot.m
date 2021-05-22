function [Ax, Fig, TitleOut, LegendOut, xlblOut, ylblOut] = s_plot(sx1, sy1, PlotProperties, Commands, Title, Legendd, SubPlot, xlbl, ylbl, FontSize, PlotType, Hidden, ScreenSize, Hold, varargin)
StrChar = @(s) isstring(s) || ischar(s);
CellFlag = false;
OneCellFlag = false;
if isempty(sx1) && isnumeric(sy1)
    if size(sy1,1)>size(sy1,2)
        sx1 = repmat((1:size(sy1,1))', size(sy1,2));
    else
        sx1 = repmat((1:size(sy1,2))', size(sy1,1)); 
    end
elseif iscell(sx1) && iscell(sy1)
    if numel(sx1)==1 && numel(sy1)>1
        sx1 = repmat(sx1, size(sy1));
    end
    if any(size(sx1)~=size(sy1))
        error('Cell Sizes do not match')
    end
    CellFlag = true;
elseif isempty(sx1) && iscell(sy1)
    CellFlag = true;
    OneCellFlag = true;
end


Fonts.title = 16;
Fonts.label = 16;
Fonts.axes = 16;
Fonts.legend = 13;
if ~isempty(FontSize)
    if iscell(FontSize)
        Fonts = update_structure(Fonts,cell2struct(FontSize(2:2:end), FontSize(1:2:end), 2));
    elseif isnumeric(FontSize) || StrChar(FontSize)
        if StrChar(FontSize), FontSize = str2double(FontSize); end
        Fonts.title = FontSize;
        Fonts.label = FontSize;
        Fonts.axes = FontSize;
        Fonts.legend = FontSize;
    elseif isstruct(FontSize)
        Fonts = update_structure(Fonts, FontSize, 'new',true);
    end
end
if isempty(PlotProperties)
    PlotProperties = {};
elseif ~iscell(PlotProperties)   
    PlotProperties = {PlotProperties};
end

if ~isa(PlotType,'function_handle')
    PlotType = @plot;
end

FigExists = false;

if isempty(SubPlot)
    if ~isempty(Hidden) && Hidden
        Fig = figure('visible','off');
    else
        Fig = figure;
    end
    SubPlot = get(Fig,'Number');
else
    FigExists = ishandle(SubPlot(1));
    if ~isempty(Hidden) && Hidden
        Fig = figure(SubPlot(1),'visible','off');
    else
        Fig = figure(SubPlot(1));
    end
    ExistingAxes = findall(Fig, 'type', 'axes');
end

if length(SubPlot)>1
    Ax(1) = subplot(SubPlot(2),SubPlot(3),SubPlot(4),'Parent',Fig);
    set(Ax(1),'FontSize',Fonts.axes);
elseif length(SubPlot)==1 && ~FigExists
    Ax(1) = axes('Parent',Fig);
    set(Ax(1),'FontSize',Fonts.axes)
elseif SubPlot(1)
    Ax(1) = ExistingAxes(end);
    if length(ExistingAxes)>=2
        Ax(2) = ExistingAxes(end-1);
    end
end
if ~isempty(Hold) && Hold
    for ax = Ax
        hold(ax, 'on')
    end
end

i = 0;
if CellFlag
    N_Ax1Plots = numel(sy1);
    for ic = 1:numel(sy1)
        if OneCellFlag
            C = sy1{ic};
            if isnumeric(C)
                if size(C,2)>size(C,1)
                    C = C.';
                end
                if isempty(PlotProperties)
                    PlotType(C(:,1), C(:,2), 'Parent',Ax(1));
                else
                    PlotType(C(:,1), C(:,2), PlotProperties{:}, 'Parent',Ax(1));
                end
            elseif isgraphics(C)
                copyobj(C,Ax(1))
                if ~isempty(C.UserData)
                    eval(['Ax(1).Children(1).' C.UserData]);
                end
            end
        else
            if isempty(PlotProperties)
                PlotType(sx1{ic}, sy1{ic}, 'Parent',Ax(1));
            else
                PlotType(sx1{ic}, sy1{ic}, PlotProperties{:}, 'Parent',Ax(1));
            end
        end
        hold(Ax(1),'on');
    end
    
    if ~isempty(varargin)
        i = 1;
        for c = varargin
            if isa(c{:},'char') || isempty(c{:})
                break
            end
            i = i+1;
        end
        if isempty(PlotProperties)
            PlotType(varargin{1:i-1},'Parent',Ax(1));
        else
            PlotType(varargin{1:i-1}, PlotProperties{:}, 'Parent',Ax(1));
        end
    end
    hold(Ax(1),'off');
else
    N_Ax1Plots = min(size(sy1));
    if ~isempty(varargin)
        i = 1;
        for c = varargin
            if isa(c{:},'char') || isempty(c{:})
                break
            end
            i = i+1;
        end
        if isempty(PlotProperties)
            PlotType(sx1, sy1, varargin{1:i-1},'Parent',Ax(1));
        else
            PlotType(sx1, sy1, varargin{1:i-1}, PlotProperties{:}, 'Parent',Ax(1));
        end
    else
        if isempty(PlotProperties)
            PlotType(sx1, sy1,'Parent',Ax(1));
        else
            PlotType(sx1, sy1, PlotProperties{:}, 'Parent',Ax(1));
        end
    end
end

br = false;
l = 1;
%% Create first axis Legend
isLegend = (iscell(Legendd) && ~isempty(Legendd{1})) || ~isempty(Legendd);
if isLegend
    if isa(Legendd,'cell')
        if size(Legendd,1)>size(Legendd,2)
            Legendd = Legendd.';
        end
        for c = Legendd
            if isa(c{:},'double')
                br = true;
                break
            end
            l = l+1;
        end
    end
    if br
        if numel(Legendd(1:l-1))<numel(findobj(Ax(1).Children, 'Type','Line', '-or', 'Type','FunctionLine'))
            Ax1Legend = Ax(1).Legend.String;
            Ax1Legend(end-N_Ax1Plots+1:end) = [];
            legend(Ax(1), [Ax1Legend Legendd(1:l-1)], 'Location','best', 'FontSize',Fonts.legend);
        else
            legend(Ax(1),Legendd{1:l-1}, 'Location','best','FontSize',Fonts.legend);
        end
    else
        if ischar(Legendd) || isstring(Legendd)
            Legendd = {Legendd};
        end
        if numel(Legendd)<numel(findobj(Ax(1).Children, 'Type','Line', '-or', 'Type','FunctionLine'))
            Ax1Legend = Ax(1).Legend.String;
            Ax1Legend(end-N_Ax1Plots+1:end) = [];
            legend(Ax(1), [Ax1Legend Legendd], 'Location','best', 'FontSize',Fonts.legend);
        else
            legend(Ax(1), Legendd, 'Location','best', 'FontSize',Fonts.legend);
        end
    end
end
if any(size(Commands))
    hold(Ax(1),'on');
    if iscell(Commands)
        for c = Commands(1,:)
            if ~isempty(c{:})
                eval(c{:});
            end
        end
    elseif isstring(Commands)
        for c = Commands(1,:)
            eval(c);
        end
    elseif ischar(Commands)
        eval(Commands);
    end
    hold(Ax(1),'on');
end

% LineStyle = {'-', '--', '-.'};
if i
    N_Ax2Plots = length(varargin(i+1:end))/2;
    if FigExists && length(ExistingAxes)>1
        Ax(2) = ExistingAxes(end-1);
    else
        Ax(2) = axes('Position',Ax(1).Position,'XAxisLocation','top','YAxisLocation','right','Color','none');
        set(Ax(2),'FontSize',Fonts.axes)
    end
	hold(Ax(2),'on');
% 	k = 0;
    for j = 0:2:(N_Ax2Plots*2 -2)
        if N_Ax2Plots==1 && min(size(varargin{i+2+j}))==1
            PlotType(varargin{i+1+j:i+2+j},'ro','Parent',Ax(2))
        else
            PlotType(varargin{i+1+j:i+2+j},'o','Parent',Ax(2))
        end
%         line(varargin{i+1+j:i+2+j},'Parent',Ax(2),'color','r','LineStyle',LineStyle{mod(k,3)+1})
%         k = k+1;
    end
        hold(Ax(2),'off');
    Ax(2).XAxis.Visible = 'off';
end

%% Create second axis Legend
if isLegend
    LegendOut(1) = Ax(1).Legend;
    if br && i
        if numel(Legendd(l+1:end))<numel(findobj(Ax(2).Children, 'Type','Line', '-or', 'Type','FunctionLine'))
            Ax2Legend = Ax(2).Legend.String;
            Ax2Legend(end-N_Ax2Plots+1:end) = [];
            legend(Ax(2), [Ax2Legend Legendd(l+1:end)], 'Location','best', 'FontSize',Fonts.legend);
        else
            legend(Ax(2),Legendd{l+1:end}, 'Location','best','FontSize',Fonts.legend);
        end
        LegendOut(2) = Ax(2).Legend;
    end
end
if ~exist('LegendOut', 'var'), LegendOut = []; end

if ~isempty(Title)
    TitleOut = title(Title,'FontSize',Fonts.title,'Parent',Ax(1));
end
if ~exist('TitleOut', 'var'), TitleOut = []; end

if isa(xlbl,'char')
    xlblOut = xlabel(xlbl,'FontSize',Fonts.label,'Parent',Ax(1));
elseif isa(xlbl,'cell')
    for i=1:length(xlbl)
        xlblOut(i) = xlabel(xlbl{i},'FontSize',Fonts.label,'Parent',Ax(i));
    end
end
if ~exist('xlblOut', 'var'), xlblOut = []; end

if isa(ylbl,'char')
    ylblOut = ylabel(ylbl,'FontSize',Fonts.label,'Parent',Ax(1));
elseif isa(ylbl,'cell')
    for i=1:length(ylbl)
        ylblOut(i) = ylabel(ylbl{i},'FontSize',Fonts.label,'Parent',Ax(i));
    end
end
if ~exist('ylblOut', 'var'), ylblOut = []; end

set(Ax(1),'FontSize',Fonts.axes);

if ~isempty(ScreenSize)
    if length(ScreenSize)==1 && ScreenSize
        Fig.WindowState = 'maximized';
    elseif length(ScreenSize)==4
        Fig.Position = ScreenSize;
    end
    
end
grid on
if any(size(Commands))
    xi = 1;
    for Axi = Ax
        hold(Axi,'on');
        if iscell(Commands)
            if size(Commands,1)>=xi
                for c = Commands(xi,:)
                    if ~isempty(c{:})
                        eval(c{:});
                    end
                end
            end
        elseif isstring(Commands)
            if size(Commands,1)>=xi
                for c = Commands(xi,:)
                    eval(c);
                end
            end
        elseif ischar(Commands)
            eval(Commands);
        end
        hold(Axi,'on');
        xi = xi+1;
    end
end
i = 1;
for ax = Ax
    ax.UserData = i;
    hold(ax, 'off')
    i = i+1;
end


end
