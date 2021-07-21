function output = absapprox(x)
beta= 1e-3;
output=sqrt(beta^2+x.^2);
% lb=-100;
% ub=100;
% np=1000;
% 
% x=linspace(lb,ub, np);
% alpha=1e3;
% beta=1e-3
% figure()
% plot(x, abs(x)); hold on 
% plot(x, 1/alpha*(log(1+exp(alpha*x))+log(1+exp(-alpha*x)))); hold on
% plot(x, (beta^2+x.^2).^0.5)
