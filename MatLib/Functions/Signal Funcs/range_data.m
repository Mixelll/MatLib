function [x,y] = range_data(x,y,varargin)
if ~isempty(varargin)
    range = varargin{:};
    logical_vec = and(x>=range(1),x<=range(2));
    x = x(logical_vec);
    y = y(logical_vec);
end
end

