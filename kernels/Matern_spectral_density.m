function S = Matern_spectral_density(f, nu, theta, D)
% Spectral density of the Matérn kernel.

rho  =  exp(theta(1));
k0  =  exp(theta(2));

% S = k0*2^D*pi^(0.5*D)*gamma(nu+D/2)*(2*nu)^nu*(2*nu/(rho^2)+4*(pi*f).^2).^(-(nu+D/2))/(gamma(nu)*rho^(2*nu));
S = k0*2^D*pi^(0.5*D)*gamma(nu+D/2)*(2*nu)^nu*(2*nu/(rho^2)+f.^2).^(-(nu+D/2))/(gamma(nu)*rho^(2*nu));

% if D == 1 && nu == 3/2
%     %S = 4*k0 * (sqrt(3)/rho)^3 * 1./((sqrt(3)/rho)^2 +(2*pi*f).^2).^2;
%     S = 4*k0 * (sqrt(3)/rho)^3 * 1./((sqrt(3)/rho)^2 + f.^2).^2; %En pulsation, pas en fréquence
% elseif D == 1  && nu == 5/2
% %    S =  k0*2*pi^(0.5)*gamma(nu+1/2)*(2*nu)^nu*(2*nu/(rho^2)+4*(pi*f).^2).^(-(nu+1/2))/(gamma(nu)*rho^(2*nu));
%    S =  k0*2*pi^(0.5)*gamma(nu+1/2)*(2*nu)^nu*(2*nu/(rho^2)+f.^2).^(-(nu+1/2))/(gamma(nu)*rho^(2*nu));
% 
% end

% S = k0*2^D*pi^(0.5*D)*gamma((3+D)/2)*(3)^(3/2)*(3/(rho^2)+4*(pi*f).^2).^(-(nu+D/2))/(gamma(nu)*rho^(2*nu));
