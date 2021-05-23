function nfig = max_fig_num()
h =  get(findobj('type','figure'),'Number');
l = length(h);
if l>1
    nfig = max([h{:}]);
elseif l==1
    nfig = h;
else
    nfig = 0;
end
end

