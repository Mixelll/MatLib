function [] = N_dim_plot(data, axis, plot, dataplot, varargin)
if ~isempty(varargin)
    datagroup = varargin{1};
else
    datagroup = [];
end
datadims = 1:ndims(data);
datadims(dataplot) = [];
datadims(datagroup) = [];
data = permute(data, [dataplot datagroup datadims]);
axis = permute(axis, [dataplot datagroup datadims]);

datasize = size(data);
dataplot = 1:length(dataplot);
datagroup = (dataplot(end) +1):(length(datagroup)+dataplot(end));

indcell = num2cell(datasize([dataplot datagroup]));
datasize([dataplot datagroup]) = [];
for i=datasize
    indcell{end+1} = ones(1,i);
end

datacells = mat2cell(data, indcell{:});
axiscells = mat2cell(axis, indcell{:});


nvc   = ndims(datacells);
vc    = ones(1, nv);
vc = num2cell(vc);
vLimc = size(datacells);
for i=1:numel(datacells)
    nv   = ndims(datacells{vc});
    v    = ones(1, nv);
    vLim = size(datacells{vc});
    
    for j=1:numel(datacells{vc})
        
        for k = 1:nv
          v{k} = v{k} + 1;
          if v{k} <= vLim(k)
            break;
          end
          v{k} = 1;
        end
    end

	for k = 1:nvc
      vc{k} = vc{k} + 1;
      if v{k} <= vLimc(k)
        break;
      end
      vc{k} = 1;
	end
end

for i=1:

    
mat2cell(data, )

sizem = size(mat3D);
figure()
hold on
C = zeros([sizem(1:2) 3]);
leg = cell(1, size(mat3D,3));
for i=1:size(mat3D,3)
    C(:,:,1) = rand();
    C(:,:,2) = rand();
    C(:,:,3) = rand();
%     C(:,:,1) = 1
%     C(:,:,2) = 0
%     C(:,:,3) = 0
    surf(mat3D(:,:,i), C);
    leg{i} = num2str(i);
    legend(leg(1:i)); 
end
hold off
end

