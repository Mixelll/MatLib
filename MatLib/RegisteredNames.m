function [StrCells, HandleCells]= RegisteredNames(varargin)
StrCells = {};
HandleCells = {};
for c = varargin
    switch c{:}
    case 'FitClasses'
        Str = {'CV_Scht_Fit_A', 'DLCP_Fit_A', 'CV_Scht_Fit_Gr1', 'CV_Scht_Fit_Gr2', 'CV_Scht_Fit_Gr3', 'CV_Scht_Fit_Gr4'};
        Handle = {@FitClasses.CV_Scht_Fit_A, @FitClasses.DLCP_Fit_A, @FitClasses.CV_Scht_Fit_Gr1, @FitClasses.CV_Scht_Fit_Gr2, @FitClasses.CV_Scht_Fit_Gr3, @FitClasses.CV_Scht_Fit_Gr4};
    end
    StrCells = [StrCells Str];
    HandleCells = [HandleCells Handle];
end
[StrCells, Indices] = unique(StrCells);
HandleCells = HandleCells(Indices);
end

