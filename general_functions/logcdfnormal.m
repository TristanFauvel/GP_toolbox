function output = logcdfnormal(x)

mask = x<=-30;

output(~mask) = log(normcdf(x(~mask)));

if any(mask(:))
    output(mask) = log(1./sqrt(2*pi))-0.5*x(mask).^2 - log(-x(mask));
end
