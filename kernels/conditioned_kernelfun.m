function [C, dC, dC_dxj] = conditioned_kernelfun(theta, kernelfun, xi, xj, training, x0, reg)
%Condition : g(x0) = 0

nj = size(xj,2);
ni = size(xi, 2);
D = size(x0,1);

if size(x0,2)>1
    warning("The function was not designed to work with more than one condition")
end
if nargout < 2
    K = kernelfun(theta, xi, xj, training, reg);
    ki = kernelfun(theta, xi, x0, training, reg);
    kj = kernelfun(theta, xj, x0, training, reg);
    k0 = kernelfun(theta, x0, x0, training, reg);
elseif nargout ==2
    [K, dK]=kernelfun(theta, xi, xj, training, reg);
    [ki, dki]=kernelfun(theta, xi, x0, training, reg);
    [kj, dkj]=kernelfun(theta, xj, x0, training, reg);
    [k0, dk0]=kernelfun(theta, x0, x0, training, reg);
    
    dC = dK - (mtimesx(dki/k0,kj,'T') - mtimesx(mtimesx(ki, dk0),kj,'T')./(k0^2)  + mtimesx(ki,dkj,'T')/k0);%derivative with respect to theta
    
elseif nargout ==3
    [K, dK, dK_dxj] =kernelfun(theta,xi, xj, training, reg);
    [ki, dki] = kernelfun(theta, x0, xi, training, reg);
    [kj, dkj, dkj_dxj] = kernelfun(theta, x0, xj, training, reg);
    [k0, dk0] =kernelfun(theta, x0, x0, training, reg);
    
    % dK_dx : ni x nj x nj x D
    % dkj_dx : 1 x nj x nj x D
%     if D>1
%         if ~isequal(size(dK_dxj), [ni,nj,nj,D]) || ~isequal(size(dkj_dxj), [1,nj,nj,D])
%             error('Error in the kernel derivatives dimensions')
%         end
%     end
    ki = ki';
    dki = permute(dki,[2 1 3]);
    kj = kj';
    dkj = permute(dkj,[2 1 3]);
    %     %%
    %     if size(xi,2)>1
    %         dkj_dx = squeeze(dkj_dx)'; % nj x D
    %         dK_dx =  squeeze(dK_dx)';
    %     else
    %                 dkj_dx = squeeze(dkj_dx)'; % nj x D
    %         dK_dx =  squeeze(dK_dx);
    %     end
    
    %%
%     dkj_dx = permute(dkj_dx,[2 1 3 4]);
    
    dC = dK - (mtimesx(dki/k0,kj,'T') - mtimesx(mtimesx(ki, dk0),kj,'T')./(k0^2)  + mtimesx(ki,dkj,'T')/k0);%derivative with respect to theta
    dC_dxj = dK_dxj - mtimesx(ki/k0,dkj_dxj);
%     if D>1
%         if ~isequal(size(dC_dxj), [ni,nj,nj,D])
%             error('Error in the kernel derivative dimensions')
%         end
%     end
    %     dC_dxj = dK_dx' - (ki/k0)*dkj_dx'; %ni x D
    
end
C = K - (ki/k0)*kj'; %ni x nj

if isequal(xi,xj)
    C = (C+C')/2; %to ensure symmetry;
end

if strcmp('reg', 'nugget')
    C = nugget_regularization(C);
end

return

