function pz = SUNpdf(z,ksi,  Omega, Omega_bar , delta, gamma, Gamma)

% gamma = s x 1;
%Gamma = s x s;
% delta= p x s;
% z = p x N;
% Omega = p x p
%ksi = p X 1 
s = numel(gamma);
p = size(z,1);
D = diag(sqrt(diag(Omega)));
m1 = delta'/Omega_bar/D; % s x p
m2 = z-ksi; % p x N

if p==1
    t1 = normpdf(z, ksi', Omega);
else
    t1 = mvnpdf(z', ksi', Omega); %argument must be N x s
end

if s== 1
    t2 = normcdf(gamma + m1*m2,0,Gamma - delta'/Omega_bar*delta);
    t3 = normcdf(gamma,zeros(1,s),Gamma);
else
%     t2 = mvnpdf((gamma + m1*m2)',0,Gamma - delta'/Omega_bar*delta); % argument must be N x s
    t2 = mvnpdf((gamma + m1*m2)',0,Gamma - delta'/Omega_bar*delta); % argument must be N x s

    t3 = mvnpdf(gamma',zeros(1,s),Gamma);    
end
t1 = t1(:);
t2= t2(:);

pz = t1.*t2/t3;