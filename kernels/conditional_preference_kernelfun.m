function [C, dC, dC_dx] = conditional_preference_kernelfun(theta, base_kernelfun, x, xp, training, reg, x0)
d =size(x,1)/2;

kernelfun = @(theta, xi, xj, training, reg) conditioned_kernelfun(theta, base_kernelfun, xi, xj, training, x0, reg);

xi= x(1:d,:);
xj = x(d+1:end,:);
xk = xp(1:d,:);
xl = xp(d+1:end,:);
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
    
    dC_dx = NaN([size(C),size(xp,2), d*2]);
    dC_dx(:,:,:,1:d) = dC_dxk; 
    dC_dx(:,:,:,(d+1):2*d) = dC_dxl;
end

if isequal(x,xp)
    C = (C+C')/2; %to ensure symmetry;
    C = nugget_regularization(C);
end
return
 