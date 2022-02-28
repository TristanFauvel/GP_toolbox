function dF = test_matrix_negderiv(f, x, eta, Icheck, verbose)
% Test derivatives when the input of the function is a matrix;

DEFAULT('eta', 1e-6);
DEFAULT('verbose', 0);

F = f(x);
dF = NaN(size(F,1), size(F,2), size(x,2), size(x,1));
count = 0;
for i = 1:size(x,2)
    for d = 1:size(x,1)
        count = count+1;
        if verbose
            display_progress(count, numel(Icheck), t0);
        end
        xdash = x;
        xdash(d,i) = x(d,i) + eta;
        Fdash = f(xdash);
        dF(:,:,i,d) = (Fdash-F)/eta;
    end
end
return
