function [s] = datasubstract(s,sub)

if size(sub,2)>1
    sb = s(:,2) - sub(:,2);
else
    sb = s(:,2) - sub;
end

sb(sb<0)=0;

s(:,2) = sb;

end

