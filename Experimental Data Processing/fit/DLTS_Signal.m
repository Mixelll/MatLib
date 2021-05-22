function [str, fun, span, FittedParamCell] = DLTS_Signal(x,y,RateWindow,N,C0,SignalMode,varargin)
isbool = @(x) x==1 || x==0;
isboolorempty = @(x) isempty(x) || x==1 || x==0;
p = inputParser;
p.KeepUnmatched=true;
p.addParameter('DecayUp', false, isboolorempty); % if C(0)<C(inf)
p.addParameter('Solve', false, isbool); % Use solve to extract tau
p.addParameter('LegendMode', 1, @isnumeric); % What to display in legend
p.addParameter('SurfaceStates', false, isbool); % Calculate surface state density instead of volume states
p.addParameter('InterfaceWidth', 0, @isnumeric); % Depletion layer width
p.parse(varargin{:});
DU = p.Results.DecayUp;
SV = p.Results.Solve;
LM = p.Results.LegendMode;
dC0 = p.Results.dC0;
SS = p.Results.SurfaceStates;
IW = p.Results.InterfaceWidth;

UserInput = @(x) ~any(strcmpi(p.UsingDefaults, x));
UU = UserInput('DecayUp') && ~isempty(DU);
if isempty(DU), DU = false; end

dvec = isnan(y);
if any(dvec)
    x = x(~dvec);
    y = y(~dvec);
end
RW1 = RateWindow(1); RW2 = RateWindow(2);
tauRW = (RW1-RW2)/log(RW1/RW2);
[~,x1Ind] = min(abs((x-RW1)));
[~,x2Ind] = min(abs((x-RW2)));
Sx1 = x(x1Ind); Sx2 = x(x2Ind);
Sy0 = y(1); Sye = y(end);
Sy1 = y(x1Ind); Sy2 = y(x2Ind);
Sy1C0 = Sy1-C0; Sy2C0 = Sy2-C0;
switch SignalMode
    case 0
        Sy1m = Sy1;
        Sy2m = Sy2;
        Signal = Sy1m - Sy2m;
        SignalModeStr = ' Absolute: C(t1)-C(t2)';
        tau = (Sx2 - Sx1)/log(abs(Sy1C0)/abs(Sy2C0));
    case 1
        Sy1m = Sy1/C0;
        Sy2m = Sy2/C0;
        Signal = Sy1m - Sy2m;
        SignalModeStr = ' Relative: C(t1)-C(t2)/C(\infty)';
        tau = (Sx2 - Sx1)/log(abs(Sy1C0)/abs(Sy2C0));
    case 2
        if dC0
            Sy1m = Sy1/dC0;
            Sy2m = Sy2/dC0;
        else
            Sy1m = Sy1/(Sy0 - Sye);
            Sy2m = Sy2/(Sy0 - Sye);
        end
        Signal = Sy1m - Sy2m;
        SignalModeStr = ' Normed: C(t1)-C(t2)/\DeltaC(0)';
        if SV
            syms tau
            f = exp(-Sx1/tau) - exp(-Sx2/tau);
            Sol = solve(f==Signal, 'real',true);
            [~,Ind] = min(abs(sol-Sx2+Sx1));
            tau = Sol(Ind);
        else
            tau = (Sx2 - Sx1)/log(Sy1m/Sy2m);
        end
    case 3
        frac = exp(-Sx1/tauRW) - exp(-Sx2/tauRW);
        Sy1m = Sy1/C0;
        Sy2m = Sy2/C0;
        if DU
            Signal = -(Sy1m - Sy2m)/frac;
            SignalModeStr = ' Trap cm^{-3}: 2*N *- C(t1)-C(t2)/C(\infty)/e^{t_1}-e^{t_2}'; 
        else
            Signal = (Sy1m - Sy2m)/frac;
            SignalModeStr = ' Trap cm^{-3}: 2*N * C(t1)-C(t2)/C(\infty)/e^{t_1}-e^{t_2}';
        end
        tau = (Sx2 - Sx1)/log(abs(Sy1C0)/abs(Sy2C0));
        
end
Nt = (Sy1 - Sy2)/(exp(-Sx1/tau)-exp(-Sx2/tau))/Sye * 2 * N;
if UU
    if DU
        Nt = -Nt;
    end
else
    SignNt = sign(Nt);
    Nt = abs(Nt);
end
if LM
    RateWindowStr = [num2str(max(x(1),RW1),1) ' to ' num2str(min(x(end),RW2),1)];
    switch LM
        case 1
            str = ['Signal= ' num2str(Signal,2)];
        case 2
            str = ['Signal= ' num2str(Signal,2) newline 'Nt= ' num2str(Nt,2) 'cm^{-3} \tau= ' num2str(tau,2) 'sec'];
        case 3
            str=['Rate Window= ' RateWindowStr  newline];
            str = [str ' Signal= ' num2str(Signal,2) SignalModeStr newline];
            str = [str 'Nt= ' num2str(Nt,2) 'cm^{-3} \tau= ' num2str(tau,2) 'sec'];
    end
else
    str = '';
end
if UU
    if DU
        fun = @(t) -Nt*exp(-t/tau)*(C0/2/N) + C0;
    else
        fun = @(t) -Nt*exp(-t/tau)*(C0/2/N) + C0;
    end
else
    if SignNt>0
        fun = @(t) Nt*exp(-t/tau)*(C0/2/N) + C0;
    else
        fun = @(t) -Nt*exp(-t/tau)*(C0/2/N) + C0;
    end
end
span = [min(x), max(x)];
FittedParamCell = {'S' ; Signal ; ['DLTS Signal' SignalModeStr]};
FittedParamCell(:, end+1) = {'Nt' ; Nt ; 'Nt - Trap Concentration'};
FittedParamCell(:, end+1) = {'tau' ; tau ; '\tau - Emission Time'};

end

