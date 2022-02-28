function dF = test_deriv(f, x, eta, Icheck, verbose)

DEFAULT('eta', 1e-6);
DEFAULT('Icheck', 1:size(x,1));
DEFAULT('verbose', 0);
 
F = f(x); 
dF = NaN(size(F,1), size(F,2), length(x));
count = 0;
for i = Icheck
    count = count+1;    
    if verbose
        display_progress(count, numel(Icheck), t0);
    end    
    xdash = x;
    xdash(i) = x(i) + eta;
    Fdash = f(xdash);    
    dF(:,:,i) = (Fdash-F)/eta;    
end 
return
