function p = tnormpdf(x, mu, sigma, a, b)

p = (1./sigma).*normpdf((x-mu)./sigma)./(normcdf((b-mu)./sigma)-normcdf((a-mu)./sigma));
