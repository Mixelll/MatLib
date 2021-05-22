function EnlargeIcons(Icons,Size)
for i = 1:length(Icons)
    if isequal(class(Icons(i)),'matlab.graphics.primitive.Line')
        Icons(i).MarkerSize = Size;
    end
end
end

