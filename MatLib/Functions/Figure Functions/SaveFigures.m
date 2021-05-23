function SaveFigures(Figures, SaveDir, CellGeneratorFuncs, SaveFormat, varargin)
CharOrString = @(s) ischar(s) || isstring(s);
p = inputParser;
p.KeepUnmatched=true;
p.addParameter('Regexp', '',CharOrString)
p.addParameter('N', 1, @isnumeric)
p.addParameter('After', false)
p.addParameter('Max', 60, @isnumeric)
p.addParameter('Append', {})
if length(varargin)==1
    p.parse
else
    p.parse(varargin{:})
end
exp = p.Results.Regexp;
N = p.Results.N;
Af = p.Results.After;
Max = p.Results.Max;
Ap = p.Results.Append;
if ~iscell(CellGeneratorFuncs)
    CellGeneratorFuncs = {CellGeneratorFuncs};
end

if ~exist(SaveDir, 'dir')
    mkdir(SaveDir)
end
Append = '';
if length(varargin)==1
    if ischar(varargin{1})
        AppendIter = varargin;
    elseif iscell(varargin{1}) || isstring(varargin{1})
        AppendIter = varargin{1};
    end
elseif iscell(Ap) || isstring(Ap)
    AppendIter = Ap;
elseif ischar(Ap)
    AppendIter = {Ap};
end
for c = AppendIter
    ci = c{:};
    if ~(isstring(ci) || ischar(ci))
        Append = [Append sprintf('%g ',x)];
    else
        Append = [Append ci];
    end
end
if length(Figures)>size(Figures,2)
    Figures = Figures';
end
Append = regexprep(Append, {'\','/','{','}'},{'','','',''});

SaveFormat = strrep(SaveFormat,'.','');
% PropertyCell = fn_struct2cell(PropertyStruct);
for f = Figures
    AppendProp = '';    
    for c = CellGeneratorFuncs
        cf = c{:};
        if isa(cf, 'function_handle')
            Prop = cf(f);
            if ~(iscell(Prop) || isstring(Prop))
                Prop = {Prop};
            end
            for p = Prop
                if iscell(p)
                    p = p{:};
                end
                if ~(isstring(p) || ischar(p))
                    p =  sprintf('%g ',p);
                end
                if isempty(exp)
                    fi1 = find(p==' ', N);
                    fi2 = find(p==',', N);
                    fi3 = find(p=='-', N);
                    fi4 = find(p==newline, N);
                    Fi = [];
                    for fic = {fi1, fi2, fi3, fi4}
                        fi = fic{:};
                        if ~isempty(fi)
                            Fi = min([Fi fi(end)]);
                        end
                    end
                else
                    [fis,fie] = regexp(p, exp);
                    if Af
                        Fi = fie(min(N,length(fie)));
                    else
                        Fi = fie(min(N,length(fis)));
                    end
                end
                if ~isempty(Fi)
                    if Af
                        p = p(Fi+1:end);
                        p = p(length(p)-min(Max,length(p))+1:end);
                    else
                        p = p(1:Fi-1);
                        p = p(1:min(Max,length(p)));
                    end
                else
                    p = p(1:min(Max,length(p)));
                end 
                AppendProp = [AppendProp ' ' p];
            end
        end
    end
    
    AppendProp = regexprep(AppendProp, {'\','/','{','}'},{'','','',''});
    SavePath = [SaveDir '\' AppendProp Append];
    for s = SaveFormat
        Pa = [SavePath '.' s{:}];
        if ~isfile(Pa)
            saveas(f, Pa)
        else
            P_i = 1;
            while isfile(Pa)
                Pa = [SavePath ' ' num2str(P_i) '.' s{:}];
                P_i = P_i+1;
            end
            saveas(f, Pa)
        end         
    end
end
end

