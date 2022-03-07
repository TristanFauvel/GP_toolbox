function [g, dgdx, dgdmu, dgdsigma] =  Gaussian_fun(x, mu, sigma)

g= normpdf(x,mu,sigma); %1/(sqrt(2*pi)*sigma)*exp(-(x-mu).^2/(2*sigma.^2))
dgdx =-(x-mu)/(sigma.^2).*g;
dgdmu= -dgdx;
dgdsigma = g.*(((x-mu)/sigma).^2-1)./sigma; %derivative of the gaussian
