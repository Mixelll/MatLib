fileprop = {'txt', '	', [0 0 100 1]}; %{'format', 'delimeter', [row_shift column_shift OPTIONAL_to_row OPTIONAL_to_column]}
path = 'D:\Cascade 301\Collective Data From Samples\gr-Si schottky P-type non-implanted b1 (big box)\right 2nd from top (center)\CV\4280 1MHz\100\try';
index = '';
saveop = 2;% 0 - don't save. 1 - save ak cell array. 2 - save an individual matrices
[CVcells, ref] = data_folder_read(path, fileprop{:}, '', index, saveop); %(path, {type, delimeter, RC}, sortby, ind, saveop)