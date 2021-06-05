classdef CV_Scht_Fit_Gr2 < matlab.mixin.SetGet
	properties
        T = 296;
        Range = [-inf inf]
        Area = (150e-4)^2 * pi
        RelativePermitivitty = 11.68
        Doping = sym('N')
        Doping_Lim = [0 1e17 inf]
        Vb = sym('Vb')
        Vb_Lim = [0 0.5 10]
        IdealityFactor = sym('n')
        Ideality_Factor_Lim = [0.01 1 100]
        Graphene_Doping = sym('n0');
        Graphene_Doping_Lim = [0 9e12 inf];
        v_fermi_Gr = 1e8;
        Graphene_Vdrop_Coeff = 'Dependent on n0';
        Graphene_Vdrop_Coeff_Lim = [-10 0.5 10 ];
        FitProperties = {'Robust','Bisquare'}
        FitPlotProperties = {}
        LegendProp = {'Location','Best'}
    end
	methods 
        function Fit(o, x,y, Target, varargin)
            DefNameValueTemp = fn_struct2cell(C_schot_fit2_B2('parser','','','','','','','',''));
            DefNameValue = [DefNameValueTemp(4,:) ; DefNameValueTemp(2,:)];
%             LimChanged = @(x,name) ~isempty([DefNameValue{2,strcmpi(DefNameValue(1,:),name)}]) && ~isequal(x,[DefNameValue{2,strcmpi(DefNameValue(1,:),name)}]);
            LimChanged = @(new,def) ~isempty(new) && ~isequal(new,def);
            Limits = {};
            FN = fieldnames(o);
            for c = DefNameValue
                if any(strcmp(c{1},FN))
                    if LimChanged(o.(c{1}), c{2}), Limits(end+1:end+2) = {c{1} o.(c{1})}; end
                end
            end
            [FitLeg, fun, span] = C_schot_fit2_B2(x,y,o.Area,o.RelativePermitivitty,o.Doping,o.Vb,o.IdealityFactor,o.Graphene_Vdrop_Coeff,o.Graphene_Doping,Limits{:});
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
            [~, OK, Pairs] = StructrureFieldsMenu(o,@parse_num_cell_sym2char,@parse_str2num_cell_sym,['Input Fit Parameters for model ' C_schot_fit2_B2('model','','','','','','','','')]);
            if OK
                set(o, Pairs{:})
            end
        end
	end
end