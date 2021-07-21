function S =SE_spectral_density(f, nu, theta, D)
% Spectral density of the Matérn kernel.

rho  =  exp(theta(1));
k0  =  exp(theta(2));

S = k0*(2*pi)^(D/2)*rho^D*exp(-rho^2*f.^2/2);

% S = k0*2^D*pi^(0.5*D)*gamma((3+D)/2)*(3)^(3/2)*(3/(rho^2)+4*(pi*f).^2).^(-(nu+D/2))/(gamma(nu)*rho^(2*nu));
