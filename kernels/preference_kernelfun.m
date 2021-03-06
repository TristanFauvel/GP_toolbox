function [C, dC, dC_dx] = preference_kernelfun(theta, kernelfun, x0, x, training, regularization)

d =size(x0,1)/2;

xi= x0(1:d,:);
xj = x0(d+1:end,:);
xk = x(1:d,:);
xl = x(d+1:end,:);
reg = 'no'; % do not regularize the base kernels
if nargout < 2
    Kik = kernelfun(theta, xi, xk, training, reg);
    Kjl = kernelfun(theta, xj, xl, training, reg);
    Kil = kernelfun(theta, xi, xl, training, reg);
    Kjk = kernelfun(theta, xj, xk, training, reg);
    C = Kik + Kjl - Kil - Kjk;
elseif nargout ==2
    [Kik, dKik] = kernelfun(theta, xi, xk, training, reg);
    [Kjl, dKjl] = kernelfun(theta, xj, xl, training, reg);
    [Kil, dKil] = kernelfun(theta, xi, xl, training, reg);
    [Kjk, dKjk] = kernelfun(theta, xj, xk, training, reg);
    C = Kik + Kjl - Kil - Kjk;
    dC = dKik + dKjl - dKil - dKjk;
elseif nargout ==3
    [Kik, dKik, dKik_dxk] = kernelfun(theta, xi, xk, training, reg);
    [Kjl, dKjl, dKjl_dxl] = kernelfun(theta, xj, xl, training, reg);
    [Kil, dKil, dKil_dxl] = kernelfun(theta, xi, xl, training, reg);
    [Kjk, dKjk, dKjk_dxk] = kernelfun(theta, xj, xk, training, reg);
    C = Kik + Kjl - Kil - Kjk;
    dC = dKik + dKjl - dKil - dKjk; 
    dC_dxk =  dKik_dxk - dKjk_dxk;
    dC_dxl = dKjl_dxl - dKil_dxl;
    
    dC_dx = NaN([size(C),size(x,2), d*2]);
    dC_dx(:,:,:,1:d) = dC_dxk; 
    dC_dx(:,:,:,(d+1):2*d) = dC_dxl;
end

if isequal(x0,x)
    C = (C+C')/2; %to ensure symmetry;
end

if strcmp(regularization, 'nugget')
C = nugget_regularization(C);
end

return
