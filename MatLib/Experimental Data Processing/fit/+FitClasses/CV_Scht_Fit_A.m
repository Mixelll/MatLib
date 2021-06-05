classdef CV_Scht_Fit_A < matlab.mixin.SetGet
	properties
        Range = [-inf inf]
        Area = (150e-4)^2 * pi
        RelativePermitivitty = 11.68
        RelativePermitivitty_Limits = [0 11.68 inf]
        Doping = sym('N')
        Doping_Limits = [0 1e17 inf]
        Vb = sym('Vb')
        Vb_Limits = [0 0.5 10]
        IdealityFactor = sym('n')
        IdealityFactor_Limits = [0.01 1 100]
        FitProperties = {}
        FitPlotProperties = {}
        LegendProp = {'Location','Best'}
    end
	methods 
        function Fit(o, x,y, Target, varargin)
            LimDef = [{'es'; [0 11.68 inf]} {'N'; [0 1e17 inf]} {'Vb'; [0 0.5 10]} {'n'; [0.01 1 100]}];
            LimChanged = @(x,name) any(x~=LimDef{2,strcmp(LimDef(1,:),name)});
            Limits = {};
            if LimChanged(o.RelativePermitivitty_Limits, 'es'), Limits(end+1:end+2) = {'es' o.RelativePermitivitty_Limits}; end
            if LimChanged(o.Doping_Limits, 'N'), Limits(end+1:end+2) = {'N' o.Doping_Limits}; end
            if LimChanged(o.Vb_Limits, 'Vb'), Limits(end+1:end+2) = {'Vb' o.Vb_Limits}; end
            if LimChanged(o.IdealityFactor_Limits, 'n'), Limits(end+1:end+2) = {'n' o.IdealityFactor_Limits}; end
            [FitLeg, fun, span] = C_schot_fit_A(x,y,o.Range,o.Area,o.RelativePermitivitty,o.Doping,o.Vb,o.IdealityFactor,o.FitProperties,Limits{:});
            hold(Target, 'on')
            if ~isempty(varargin)
                legend(Target, varargin{:});
            elseif isempty(isempty(Target.Legend))
                legend(Target, 'Measured Data');
            end
            Lines = findobj(Target.Children, 'Type','Line');
            
            if ~isempty(Lines)
                fplot(fun, span, o.FitPlotProperties{:}, 'Parent',Target, 'Color',Lines(1).Color);
            else
                fplot(fun, span, o.FitPlotProperties{:}, 'Parent',Target);
            end
            legend(Target,[Target.Legend.String(1:end-1) {[FitLeg newline 'Model: ' Target.Legend.String{end}]}])
            if ~isempty(o.LegendProp)
                legend(Target, o.LegendProp{:});
            end
            hold(Target, 'off')
        end
        function OK = Menu(o)
            [~, OK, Pairs] = StructrureFieldsMenu(o,@parse_num_cell_sym2char,@parse_str2num_cell_sym,['Input Fit Parameters for model ' C_schot_fit_A('model','','','','','','','','')]);
            if OK
                set(o, Pairs{:})
            end
        end
	end
end