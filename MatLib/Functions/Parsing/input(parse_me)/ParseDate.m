function [TimeStamp, DateOut] = ParseDate(Date)
if ~iscell(Date)
    Date = {Date};
end
TimeStamp = zeros(1, length(Date));
DateOut = NaT(1, length(Date));
D_i = 1;
for Dc = Date
    D = Dc{:};
M = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
m = {'jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'};
Month = 0;
for i=1:12
    if any(strfind(D, M{i})) || any(strfind(D, m{i}))
        Month = i;
        break
    end
end
TimeInd = regexp(D, ':');
TimeSeconds = str2double(D(TimeInd(1)-2:TimeInd(1)-1))*3600 + str2double(D(TimeInd(2)-2:TimeInd(2)-1))*60 + str2double(D(TimeInd(2)+1:TimeInd(2)+2));
D(TimeInd(1)-2:TimeInd(2)+2) = [];

YearInd = regexp(D,'\d{4}');
Year = str2double(D(YearInd:YearInd+3));
D(YearInd:YearInd+3) = [];

DayOfMonthInd = regexp(D,'\d{1,2}');
DayOfMonth = str2double(D(DayOfMonthInd:DayOfMonthInd+1));
if isnan(DayOfMonth)
    DayOfMonth = str2double(D(DayOfMonthInd));
end
YMD_UNIX = posixtime(datetime([Year Month DayOfMonth]));
TimeStamp(D_i) = YMD_UNIX + TimeSeconds;
DateOut(D_i) = datetime(TimeStamp, 'ConvertFrom', 'posixtime');
D_i = i+1;
end
end

