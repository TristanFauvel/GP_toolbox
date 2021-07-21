function [phi, dphidx] = quadratic_features(x)

N=size(x,2);
Nd = size(x,1)+1;
Nfeatures = (Nd+1)*Nd/2;

phi = zeros(Nfeatures,N);

id= tril(ones(Nd,Nd))-eye(Nd);

for i = 1:N
    m =triu([x(:,i);1]*[x(:,i);1]');
    m(logical(id(:)))=[];
    phi(:,i) = m;
end


if nargout>1
    diag_id = eye(Nd);
    dphidx = zeros(Nfeatures, N,Nd-1);
    for i = 1:N
        for j = 1:(Nd-1)
            m= zeros(Nd);
            m(1:Nd-1,j)= x(:,i);
            m(j,1:Nd-1) = x(:,i);  
            m(j,Nd) = 1;
            m(j,j)=2*m(j,j);            
            m= triu(m);           
            m(logical(id(:)))=[];
            dphidx(:,i,j) = m;
        end
    end    
end


% if nargout>1
%     diag_id = eye(Nd);
%     dphidx = zeros(Nfeatures, Nd-1, N);
%     for i = 1:N
%         for j = 1:(Nd-1)
%             m =triu(([x(:,i);1]*ones(1,Nd))');
% %             m =triu((ones(Nd,1))*[x(:,i);1]');
%             m(logical(eye(Nd))) = m(logical(eye(Nd)))*2;
%             m(setdiff(1:Nd,j),setdiff(1:Nd,j))= 0;
%             m(logical(id(:)))=[];
%             dphidx(:,j,i) = m;
%         end
%     end    
% end

% if nargout>1
%     diag_id = eye(Nd);
%     dphidx = zeros(Nfeatures, Nd-1, N);
%     for i = 1:N       
%         m =triu(([x(:,i);1]*ones(1,Nd))');
%         m(end,end)=0; %derivative of the bias
%         m(logical(id(:)))=[];
%         dphidx(:,i) = m;
%     end
%     
% end
