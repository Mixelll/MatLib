function save_plot_n(fig, PathFileName, SaveFormat)

if ~exist(fileparts(PathFileName), 'dir')
    mkdir(fileparts(PathFileName))
end
SaveFormat = strrep(SaveFormat,'.','');

j=1;
for c = SaveFormat
    i=1;
    while true
      if ~isfile([PathFileName '_' num2str(i) '.' c{:}])
          break
      else
          i=i+1;
      end
    end
    j = max(j,i);
end

for c = SaveFormat
    saveas(fig, [PathFileName '_' num2str(j) '.' c{:}]);
end
end

