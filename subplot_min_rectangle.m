function [subrow,subcol] = subplot_min_rectangle(ns)
subvec = 1:10; % max 10 by 10 graphs
submat = subvec'*subvec; % create the matrix
subrestrict = 2; % disallow long/narrow subplot schemes (only keep subrestrict diagonals around the central diagonal)  
submat(logical(tril(submat,-subrestrict))) = 0; % set disallowed diagonals to zero 
submat(logical(triu(submat,subrestrict))) = 0; % set disallowed diagonals to zero 
[subdiv, ~] = find(submat==min(submat((submat - ns)>=0))); % find size of smallest area allowed rectable that can fit the graphs
subdiv = sort(subdiv); % sort to enable to rectangular to be longer on the horizontal axis (like computer screens)
subrow = subdiv(floor((length(subdiv)+1)/2)); % generate row size
subcol = subdiv(ceil((length(subdiv)+1)/2)); % generate column size
end

