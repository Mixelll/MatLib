function out = LogSpace(a,b,n)
% if isempty(varargin)
%     Left2Right = true;
% else
%     Left2Right = varargin{1};
% end
% if Left2Right
%     out = log10(logspace(a,b,n));
% else
out = logspace(log10(a),log10(b),n);
end

