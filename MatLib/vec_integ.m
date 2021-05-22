function syinteg = vec_integ(sx, sy)

if size(sx, 2)==1
    sx = repmat(sx, 1, size(sy, 2));
end

lsy = size(sy,1);
syinteg = zeros(lsy-1, size(sy,2));
diffx = sx(2:end,:) - sx(1:end-1,:);
avgy = (sy(2:end,:) + sy(1:end-1,:))/2;
syinteg(1,:) = avgy(1,:).*diffx(1,:);

for i=2:lsy-1
    syinteg(i,:) = syinteg(i-1,:) + avgy(i,:).*diffx(i,:);
end
    
end

