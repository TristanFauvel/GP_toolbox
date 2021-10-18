function S = Matern_spectral_density(f, nu, theta, D)
% Spectral density of the Matérn kernel.

rho  =  exp(theta(1));
k0  =  exp(theta(2));

S = k0*2^D*pi^(0.5*D)*gamma(nu+D/2)*(2*nu)^nu*(2*nu/(rho^2)+f.^2).^(-(nu+D/2))/(gamma(nu)*rho^(2*nu));

 