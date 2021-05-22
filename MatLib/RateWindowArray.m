function RateWindows = RateWindowArray(Ratio, t1Start, t2End, N_Windows, DensityFunction)
if isempty(DensityFunction)
    DensityFunction = @linspace;
end
t1Array = DensityFunction(t1Start, t2End/Ratio, N_Windows);
RateWindows = [t1Array; Ratio*t1Array];
end

