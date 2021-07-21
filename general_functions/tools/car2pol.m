function [r, theta] = car2pol(x)

r= sqrt(sum(x.^2));
nd= size(x,1);
N = size(x,2);
theta = NaN(nd-1,N);

for i=1:nd-2
    theta(i,:) = acos(x(i,:)./sqrt(sum(x(i:nd,:).^2)));
end
i=nd-1;
for k =1:N
    if x(nd,k)<0
        theta(i,k) = 2*pi - acos(x(i,k)./sqrt(sum(x(i:nd,k).^2)));        
    else
        theta(i,k) = acos(x(i,k)./sqrt(sum(x(i:nd,k).^2)));
    end
end
return

