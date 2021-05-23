function varargout = demoivre(z,n,varargin)


r=abs(z);
th=angle(z);

rn = r^(1/n);
angs = (th + 2*(0:n-1)*pi)/n;
zn = rn .* (cos(angs) + 1i*sin(angs));

varargout  = num2cell(zn);
if ~isempty(varargin)
    varargout = varargout(varargin{:}+1);
end


    
end

