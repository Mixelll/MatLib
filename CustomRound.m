function Rounded = CustomRound(Num,Bra,varargin) % Number to round (arrays supported), round to (arrays supported), round up or down 
Size = size(Num);
ByElementFlag = false;
CharOrString = @(s) ischar(s) || isstring(s);
if numel(Bra)==1
    Bra = repmat(Bra, size(Num));
end
if ~sign(Num)==0
    Bra = sign(Num).*abs(Bra);
end

if ~isempty(varargin)
    f = varargin{:};
    if isequal(f,@ceil)
        Up = true;
    elseif isequal(f,@floor)
        Up = false;
    elseif CharOrString(f) && (strcmpi(f,'up') || strcmpi(f,'ceil'))
        Up = true;
    elseif CharOrString(f) && (strcmpi(f,'down') || strcmpi(f,'floor'))
        Up = false;
    else
        ByElementFlag = true;
    end
else
    ByElementFlag = true;
end

if ByElementFlag
    Up = false(1, numel(Num));
    for i = 1:numel(Num)
        Mod = abs(mod(Num(i),Bra(i)));
        if Mod<abs(Bra(i))/2
            Up(i) = false;
        else
            Up(i) = true;
        end
    end
else
    Up = repmat(Up,1,numel(Num));
end
Rounded = nan(1, numel(Num));
for i = 1:numel(Num)
    if Up(i)   
        Rounded(i) = ceil(Num(i)/Bra(i))*Bra(i); % Up
    else
        Rounded(i) = floor(Num(i)/Bra(i))*Bra(i); % Down
    end
end
Rounded = reshape(Rounded, Size);
% if strcmp(Rou, 'up')    
%     Rounded = floor(Num) + ceil( (Num-floor(Num))/Bra) * Bra; % Up
% else
%     Rounded = floor(Num) + floor( (Num-floor(Num))/Bra) * Bra; % Down
% end
end

