%% Same thing with nd = 2
add_gp_module;
figure_path = '/home/tfauvel/Documents/PhD/Figures/Thesis_figures/Chapter_1/';


rng(12) 
n=50;

nd = 2;
x = linspace(0,5,n);

[p,q] = meshgrid(x,x);
x2d = [p(:), q(:)]';

kernelfun = @ARD_kernelfun;
kernelname = 'ARD';
theta.cov= [0; 0; 0];
gen_kernelfun = @ARD_kernelfun;

g =  mvnrnd(constant_mean(x2d,0), gen_kernelfun(theta.cov(1:3), x2d,x2d)); %generate a function
g=g-mean(g);

fig = figure();
fig.Color = [0,0,0];
s = surf(reshape(-g, n,n));
s.EdgeColor = 'none';
set(gca,'visible','off')

figname  = 'GP_2d';
folder = [figure_path,figname];
savefig(fig, [folder,'/', figname, '.fig'])
exportgraphics(fig, [folder,'/' , figname, '.pdf']);
exportgraphics(fig, [folder,'/' , figname, '.png'], 'Resolution', 300);


sigma = sqrt(exp(theta.cov(end)));

y = g + sigma*randn(1,n^nd);

colo= othercolor('GnBu7');
Fontsize =14;
   
N=50;
idxtrain = randsample(size(x2d, 2),N);
xtrain = x2d(:,idxtrain);
ytrain = y(idxtrain)';

theta.mean =0;
[posterior_mean, posterior_variance, ~, ~, Sigma2_y]=prediction(theta, xtrain, ytrain', xtrain, kernelfun, @constant_mean);

m = 10;
sample_g_x2d = NaN(n^nd, m);
for i =1:m
   sample= sample_GP(theta.cov, xtrain, ytrain, kernelname,approximation_method, decoupled_bases, nfeatures, kernelfun);
 sample_g_x2d(:,i) = sample(x2d);
end

[posterior_mean, posterior_variance]=prediction(theta, xtrain, ytrain', x2d, kernelfun, @constant_mean);

h=figure(2);
h.Color =  [1 1 1];
h.Name = 'Value  function';
subplot(1,2,1)
imagesc(x, x, reshape(posterior_mean,n,n))
xlabel('x','Fontsize',Fontsize)
ylabel('x','Fontsize',Fontsize)
set(gca,'YDir','normal')
pbaspect([1 1 1])
title('True posterior mean')
colorbar
subplot(1,2,2)
imagesc(x, x, reshape(mean(sample_g_x2d, 2),n,n))
xlabel('x','Fontsize',Fontsize)
ylabel('x','Fontsize',Fontsize)
set(gca,'YDir','normal')
pbaspect([1 1 1])
title('Sample posterior mean')
colorbar()
    
h=figure(3);
h.Color =  [1 1 1];
% h.Name = 'Preference function';
subplot(1,2,1)
imagesc(x, x, reshape(posterior_variance,n,n))
xlabel('x','Fontsize',Fontsize)
ylabel('x','Fontsize',Fontsize)
set(gca,'YDir','normal')
pbaspect([1 1 1])
title('True posterior variance')
colorbar()

subplot(1,2,2)
imagesc(x, x, reshape(var(sample_g_x2d'),n,n))
xlabel('x','Fontsize',Fontsize)
ylabel('x','Fontsize',Fontsize)
set(gca,'YDir','normal')
pbaspect([1 1 1])
title('Sample posterior variance')
colorbar()

