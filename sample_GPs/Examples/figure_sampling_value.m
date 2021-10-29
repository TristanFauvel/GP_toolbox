clear all
figure_path = '/home/tfauvel/Documents/PhD/Figures/Thesis_figures/Chapter_1/';
graphics_style_paper
close all
rng(1)

letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

link = @normcdf; %inverse link function

%% Define the range of parameters
n = 30;
lb = 0; ub =1;
x = linspace(lb, ub, n);
d =1;
ntr = 5;

[p,q]= meshgrid(x);
x2d = [p(:), q(:)]';

x0 = x(:,1);

modeltype = 'exp_prop'; % Approximation method
base_kernelfun =  @Matern52_kernelfun;%kernel used within the preference learning kernel, for subject = computer
kernelname = 'Matern52';
approximation.method = 'RRGP'; %%RRGP
approximation.nfeatures = 256;

theta.cov = [log(1/10),0];
theta.mean = 0;

approximation.decoupled_bases = 1;
condition.x0 = 0;
condition.y0 = 0;

kernelfun = @(theta, xi, xj, training, regularization) conditional_preference_kernelfun(theta, base_kernelfun, xi, xj, training,regularization, condition.x0);
link = @normcdf; %inverse link function for the classification model


regularization = 'nugget';
meanfun = 0;
type = 'preference';
hyps.ncov_hyp =2; % number of hyperparameters for the covariance function
hyps.nmean_hyp =0; % number of hyperparameters for the mean function
hyps.hyp_lb = -10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
hyps.hyp_ub = 10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
D = 1;
 
model = gp_preference_model(D, meanfun, kernelfun, regularization, hyps, lb, ub, type, link, modeltype, kernelname,  condition, base_kernelfun);




g = mvnrnd(zeros(1,n),base_kernelfun(theta.cov, x, x, 'false', regularization));
g = g-g(1);

figure();
plot(g);

f = g'-g;
f= f(:);

nsamp= 5;
rd_idx = randsample(size(x2d,2), nsamp, 'true');
xtrain= x2d(:,rd_idx);
ytrain= f(rd_idx);
ctrain = link(ytrain)>rand(nsamp,1);



[~,  mu_g, sigma2_g, Sigma2_g, dmuc_dx, dmuy_dx, dsigma2y_dx, dSigma2y_dx, var_muc, dvar_muc_dx,post] = model.prediction(theta, xtrain(:,1:ntr), ctrain(1:ntr), [x; x0*ones(1,n^d)], []);
mu_g = -mu_g; %(because prediction_bin considers P(x1 > x2);



%% Find the true global optimum of g
[gmax, id_xmax] = max(g);
xmax = x(id_xmax);

nsamps = 1000;
samples_g = NaN(nsamps, n);
samples_prior = NaN(nsamps, n);
samples_f =  NaN(nsamps, n*n);
updates = NaN(nsamps, n);
 %%
for j = 1:nsamps
    [sample_f, samples_g(j,:), decomposition] = sample_preference_GP(x, theta, xtrain, ctrain, model, approximation, post);
    samples_prior(j,:) = decomposition.sample_prior(x);
    updates(j,:) = decomposition.update_2([x;x0.*ones(D,size(x,2))]);
    samples_f(j,:) = sample_f;
end


legend_pos = [-0.14,1.0];

mr =1;
mc = 3;
fig=figure('units','centimeters','outerposition',1+[0 0 16 height(mr)]);
fig.Color =  background_color;
i = 0;
tiledlayout(mr,mc, 'TileSpacing' , 'tight', 'Padding', 'tight')

nexttile();
i=i+1;
plot_gp(x, mu_g, sigma2_g, C(1,:), linewidth); 
xlabel('$x$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1])
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)

text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
box off
xlabel('$x$')

nexttile();
i=i+1;
plot(x, samples_prior(j,:), 'Color',  C(1,:),'LineWidth', linewidth); hold on ;
plot(x, updates(j,:), 'Color', C(2,:),'LineWidth', linewidth); hold on ;
plot(x, samples_g(j,:), 'Color', C(3,:),'LineWidth', linewidth); hold off ;
legend('Prior sample', 'Update', 'Posterior sample')
box off;
legend boxoff
xlabel('$x$')

text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
set(gca, 'Xlim', [0,1], 'Xtick', [0,0.5,1],'Fontsize', Fontsize);

nexttile();
i=i+1;
imagesc(x,x,reshape(samples_f(j,:),n,n))
pbaspect([1,1,1])
ylabel('$x''$')
xlabel('$x$')
colormap(cmap)
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
set(gca,'YDir','normal', 'Fontsize', Fontsize)



figure()
crosscov= cov(samples_g);

cl= [min([crosscov(:);Sigma2_g(:)]), max([crosscov(:);Sigma2_g(:)])];

subplot(1,2,1)
imagesc(x,x,Sigma2_g, cl)
pbaspect([1,1,1])
title('Posterior covariance')
set(gca,'YDir','normal')
ylabel('$x''$')
xlabel('$x$')
colormap(cmap)

subplot(1,2,2)
imagesc(x,x,crosscov)
pbaspect([1,1,1])
title('Sample covariance')
set(gca,'YDir','normal')
ylabel('$x''$')
xlabel('$x$')

colormap(cmap)



figure();
plot(samples_prior');


figname  = 'value_sampling_GP';
folder = [figure_path,figname];
savefig(fig, [folder,'/', figname, '.fig'])
exportgraphics(fig, [folder,'/' , figname, '.pdf']);
exportgraphics(fig, [folder,'/' , figname, '.png'], 'Resolution', 300);

