function S = SE_ARD_spectral_density(f, theta, D)
% Spectral density of the SE ARD kernel.
nd = numel(theta)-1;
lambda       =  exp(theta(1:nd));
lambda = lambda(:)';
k0           =  exp(theta(nd+1));

S = k0*sqrt(2*pi)^D./sqrt(prod(lambda)).*exp(-0.5*sum(f.^2./lambda(:)',2));
