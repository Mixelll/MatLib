function [] = surf3d(in1)
mat3D = in1;

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

