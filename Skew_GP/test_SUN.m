upath = userpath;
addpath(genpath([upath,'/BO_toolbox']))
addpath(genpath([upath,'/GP_toolbox']))
addpath(genpath([upath,'/Preference_Based_BO_toolbox']))

% gamma = 1 x s;
%Gamma = s x s;
% delta= p x s;
% z = p x N;
% Omega = p x p
%ksi = p X 1 (ou 1X p?)

%% Fig 1. a1
clear all

ksi = 0;
Omega_bar = 1;
D = 1;
Omega = D*Omega_bar*D;

gamma= 0;
Gamma = 1;
s =1;

z = linspace(-3,3,100);
delta_range = [0.8, -0.8,0.3,-0.3];
figure()
legends = {};
for i = 1:4
    delta =delta_range(i);
    pz = SUNpdf(z, ksi,  Omega, Omega_bar, delta, gamma, Gamma);
    plot(z, pz); hold on;
    legends{i} = ['SUN delta = ', num2str(delta)];
end
plot(z, normpdf(z));
hold off;
legends{i+1} = 'N(0,1)';
legend(legends)

%%
%% Fig 1. a2
clear all

s= 2;
p=1;
z = linspace(-3,3,100);
ksi = 0;

Gamma = [1,0.8;0.8,1];

Omega_bar = 1;
D = 1;
Omega = D*Omega_bar*D;

gamma = zeros(s,1);

figure()
legends = {};

% pz = zeros(1,100);
delta_range = [0.8,0.8,-0.8,-0.8;0.3,-0.3,-0.3,0.3];
for i = 1:4
    delta =delta_range(:,i)'; % p x s
%     for j =1:100
%         pz(j) = SUNpdf(z(j),ksi, Omega, delta, gamma, Gamma);
%     end
    pz = SUNpdf(z,ksi,  Omega, Omega_bar, delta, gamma, Gamma);
    legends{i} = ['SUN delta = ', num2str(delta(:)')];
    plot(z, pz); hold on;
end
plot(z, normpdf(z));hold off;
hold off;
legends{i+1} = 'N(0,1)';
legend(legends)


%% Fig 2. a1
clear all

s=  1;
n=1000;
zrange = linspace(-2,2,n)';
[p,q] = ndgrid(zrange);
z= [p(:),q(:)]';
p=2;

delta = [0.8;0.3];
Gamma = 1;

ksi = [0;0];

Omega_bar = [1, 0.8 ; 0.8,1];
D = eye(s);
Omega = D*Omega_bar*D;

gamma = zeros(s,1);

pz = SUNpdf(z,ksi, Omega, Omega_bar, delta, gamma, Gamma);

h = figure();
imagesc(zrange, zrange, reshape(pz, n,n))
pbaspect([1 1 1])
set(gca,'YDir','normal')
colorbar()

%% Fig 2. a2
clear all

s=  1;
n=1000;
zrange = linspace(-2,2,n)';
[p,q] = ndgrid(zrange);
z= [p(:),q(:)]';
p=2;

delta = [0.3;0.8];
Gamma = 1;

ksi = [0;0];

Omega_bar = [1, 0.8 ; 0.8,1];
D = eye(s);
Omega = D*Omega_bar*D;

gamma = zeros(s,1);

pz = SUNpdf(z,ksi, Omega, Omega_bar, delta, gamma, Gamma);

h = figure();
imagesc(zrange, zrange, reshape(pz, n,n))
pbaspect([1 1 1])
set(gca,'YDir','normal')
colorbar()

%% Fig 2. a3
clear all
s= 2;
n=1000;
zrange = linspace(-2,2,n)';
[p,q] = ndgrid(zrange);
z= [p(:),q(:)]';
p=2;

delta = [0.3, 0.5 ; 0.8, 0.9];
Gamma = [1,0.9;0.9,1];

ksi = [0;0];%ksi = p X 1 


Omega_bar = [1, 0.8 ; 0.8,1];
D = eye(s);
Omega = D*Omega_bar*D;

gamma = zeros(s,1);

pz = SUNpdf(z,ksi, Omega, Omega_bar, delta, gamma, Gamma);

figure()
imagesc(zrange, zrange, reshape(pz, n,n))
pbaspect([1 1 1])

%% Fig 2. a4
clear all

s= 2;
n=1000;
zrange = linspace(-2,2,n)';
[p,q] = ndgrid(zrange);
z= [p(:),q(:)]';
p=2;

delta = [0.3, 0.8 ; 0.8, 0.3];
Gamma = [1,-0.3;-0.3,1];

ksi = [0;0];

Omega_bar = [1, 0.8 ; 0.8,1];
D = eye(s);
Omega = D*Omega_bar*D;

gamma = zeros(s,1);

pz = SUNpdf(z,ksi, Omega, Omega_bar, delta, gamma, Gamma);

figure()
imagesc(zrange, zrange, reshape(pz, n,n))
pbaspect([1 1 1])

M= NaN(s+p);
M(1:s,1:s) = Gamma;
M(s+1:end,s+1:end) = Omega_bar; 
M(1:s,s+1:end) = delta'; 
M(s+1:end, 1:s) =  delta; 
chol(M)