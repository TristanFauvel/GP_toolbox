function x = pol2car(theta,r)
%Compute the cartesian coordinates corresponding to hyperspherical coordinates
N= size(theta,2);
nd = size(theta,1)+1;

x=NaN(nd, N);
sin_theta = sin(theta);
cos_theta = cos(theta);

% x(1,:) = r*cos_theta(1,:);

% if nd == 2
%     x(1,:)= r*cos_theta(1,:);
%     
%     x(2,:) = r*sin_theta(1,:);
% else 
%     for i = 2:nd
%         x(i,:) = r*prod(sin_theta(1:2:i,:)).*prod(cos_theta(2:2:i,:));
%     end
% end

for i = 1:nd-1
    x(i,:) = r*prod(sin_theta(1:(i-1),:)).*cos_theta(i,:);
end
x(nd,:) = r*prod(sin_theta);

return