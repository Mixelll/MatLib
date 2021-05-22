function [Indices, GroupsMean] = Groups(G,D) % Elements to group, max distance of outer elements
Size = size(G);
[Gs, Gi]= sort(G(:).');
LGs = length(Gs);
Gdi = find(diff(Gs)>D);
GroupsMean = nan(1,length(Gdi));
Gdi = [0 Gdi LGs];
Gsi = nan(1,LGs);
for i = 1:length(Gdi)-1
    GroupsMean(i) = mean(Gs(Gdi(i)+1:Gdi(i+1)));
    Gsi(Gdi(i)+1:Gdi(i+1)) = i;
end
Indices(Gi) = Gsi;
Indices = reshape(Indices, Size);
end

