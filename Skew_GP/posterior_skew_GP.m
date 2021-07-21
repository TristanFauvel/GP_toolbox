% function [mu_c, mu_y,, sigma2_y] =  prediction_skew_GP
kernelfun = @Gaussian_kernelfun;
n=100;
x= linspace(0,1,n);
theta= [5,3];
K = Gaussian_kernelfun(theta,x,x);
y = mvnrnd(zeros(1,n),K);
p= normcdf(2*y-1);
% 
% figure()
% plot(p)
% 
ndata = 50;
id_data = randsample(n,ndata);
xdata = x(id_data);
cdata = p(id_data)>rand(1,ndata);


W = diag(2*cdata-1);

%% Posterior: Theorem 1 in Benavoli et al 2020
ksi = zeros(ndata,1);
Omega = Gaussian_kernelfun(theta,xdata,xdata);
D = diag(sqrt(diag(Omega)));
Omega_bar=inv(D)*Omega*inv(D);
delta_tilde = Omega_bar*D*W';
gamma_tilde = W*ksi;
Gamma_tilde = W*Omega*W' +  eye(ndata);


%% Predictive distribution
ntest= 1;
mu_c= NaN(1,n);
mu_c_mvn= NaN(1,n);
mu_c_newmvncdf= NaN(1,n);
for id_test= 1:n
    xtest = x(id_test);
    
    ksi = zeros(ndata+ntest,1);
    Omega = Gaussian_kernelfun(theta,[xdata,xtest],[xdata,xtest]);
    D = diag(sqrt(diag(Omega)));
    Omega_bar=inv(D)*Omega*inv(D);
    W = diag([2*cdata-1, 0.5*ones(1,ntest)]);
    
    % delta_tilde = Omega_bar*D*W';
    gamma_tilde_s = W*ksi;
    Gamma_tilde_s = W*Omega*W' +  eye(ndata+ntest);
    
    mu_c_newmvncdf(id_test) = newmvncdf(-inf(ndata+ntest,1),gamma_tilde_s',Gamma_tilde_s,1000).prob/newmvncdf(-inf(ndata,1),gamma_tilde',Gamma_tilde,1000).prob;
   % mu_c_mvn(id_test) = mvncdf(gamma_tilde_s', 0, Gamma_tilde_s)/mvncdf(gamma_tilde', 0, Gamma_tilde);
%     mu_c(id_test) =mvnxpb(Gamma_tilde_s, -inf(ndata+ntest,1),gamma_tilde_s)/mvnxpb(Gamma_tilde, -inf(ndata,1), gamma_tilde);
end

figure();
plot(x, p); hold on;
plot(x, mu_c_mvn); hold on;
plot(x, mu_c_newmvncdf); hold on;
plot(x, mu_c); hold on;
scatter(xdata, cdata, 25,'k', 'filled'); hold off
legend('true', 'new mvncdf', 'matlab', 'Genz', 'data')

mu_c = mu_c_newmvncdf;
%% Comparison with Laplace approximation
mu_c_Laplace =  prediction_bin(theta, xdata, cdata, x, kernelfun, 'modeltype', 'laplace');

%% Comparison with EP
mu_c_EP =  prediction_bin(theta, xdata, cdata, x, kernelfun, 'modeltype', 'exp_prop');


figure();
plot(x, p); hold on;
plot(x, mu_c); hold on;
plot(x, mu_c_Laplace); hold on;
plot(x, mu_c_EP); hold on;
scatter(xdata, cdata, 25,'k', 'filled'); hold off
legend('True function','Skew GP', 'Laplace', 'EP', 'Data')
% 
% 
% delta = Omega*Domega*W';
% 
% gamma_tilde = W*ksi;
% 
% Gamma_tilde = NaN(ntrain+s);
% Gamma_tilde(1:ntrain,1:ntrain) = Gamma;
% Gamma_tilde((ntrain+1):end, (ntrain+1):end) = W*Omega*W' +  eye(p);
% Gamma_tilde(1:ntrain, (ntrain+1):end) = delta'*D*W';
% Gamma_tilde((ntrain+1):end, 1:ntrain) = W*D*delta;
% 
% Gamma_s = NaN(ntrain+s +ntest);
% Gamma_s(1:ntrain+s,1:ntrain+s) = Gamma_tilde;
% Gamma_s(ntrain+s+1:end
% 
% mu_y = 
% 
% mu_c = normcdf(gamma_s, Gamma_s)./normcdf(gamma_tilde, Gamma_tilde);