%% Sample GP

clear all
close all
add_gp_module;
figure_path = '/home/tfauvel/Documents/PhD/Figures/Thesis_figures/Chapter_1/';

rng(13);%12
n=1000; %resolution
x = linspace(0,1,n);
ntest=1000;
xtest = linspace(0,1,ntest);

% gen_kernelfun = @Matern52_kernelfun;
% kernelfun = @Matern52_kernelfun;
% kernelname = 'Matern52';
gen_kernelfun = @ARD_kernelfun;
kernelfun = @ARD_kernelfun;
kernelname = 'ARD';

approximation.method = 'RRGP';
link = @normcdf;

modeltype = 'exp_prop';
post = [];
regularization = 'nugget';

model.regularization = regularization;
model.kernelfun = kernelfun;
model.link = link;
model.modeltype = modeltype;
model.kernelname = kernelname;
model.D = 1;

theta.cov = [4,-1];

% theta.cov = [log(1/10),0];
theta_gen.cov = theta.cov;
theta.mean= 0;

Sigma =gen_kernelfun(theta_gen.cov, x,x, true, 'false');



% figure(); plot(link(g))


graphics_style_paper;

N=5;
num_tr = [2,4,8,16,32,64,128,256,512]; %,1024];
nxtr = numel(num_tr);
hyp = theta.cov;
m=5000;
fx = NaN(m, n);

fx_decoupled_Hilbert= NaN(m, n);
fx_weight_Hilbert= NaN(m, n);
fx_decoupled_Fourier= NaN(m, n);
nrepets = 32;
W_decoupled_Hilbert = NaN(nrepets, nxtr);
W_weight_Hilbert = NaN(nrepets, nxtr);
W_decoupled_Fourier = NaN(nrepets, nxtr);

for k =1:nrepets
    g =  mvnrnd(constant_mean(x,0), gen_kernelfun(theta_gen.cov, x,x, true, 'false')); %generate a function
    for j = 1:nxtr
        disp(['Repetition :',num2str(k), ', ', num2str(j)])
        ntr = num_tr(j);
        
        idxtrain = randsample(n,ntr, 'true');
        N = numel(idxtrain);
        xtrain = x(idxtrain);
        y_data = g(idxtrain)';
        ctrain = link(y_data)>rand(N,1);
        
        post = prediction_bin(hyp, xtrain, ctrain, [], model, post);
        
        [mu_c,  mu_f, sigma2_f, Sigma2_f] = prediction_bin(hyp, xtrain, ctrain, xtest, model, post);
        
        %%
        rng(k)
        approximation.decoupled_bases = 1;
        approximation.method = 'RRGP';
        nfeatures = 32;
        approximation.nfeatures = nfeatures;
        [approximation.phi, approximation.dphi_dx] = sample_features_GP(theta.cov, model, approximation);
        for i =1:m
            gs = sample_binary_GP_precomputed_features(xtrain, ctrain, theta.cov, model, approximation, post);
            fx_decoupled_Hilbert(i, :)=gs(xtest);
        end
        W_decoupled_Hilbert(k,j) = Wasserstein2(mu_f, Sigma2_f, fx_decoupled_Hilbert);
        
        %%
        rng(k)
        approximation.decoupled_bases = 0;
        approximation.method = 'RRGP';
        nfeatures = ntr + 32;
        approximation.nfeatures = nfeatures;
        [approximation.phi, approximation.dphi_dx] = sample_features_GP(theta.cov, model, approximation);
        for i =1:m
            gs = sample_binary_GP_precomputed_features(xtrain, ctrain, theta.cov, model, approximation, post);
            fx_weight_Hilbert(i, :)=gs(xtest);
        end
        W_weight_Hilbert(k,j) = Wasserstein2(mu_f, Sigma2_f, fx_weight_Hilbert);
        %%
        rng(k)
        approximation.decoupled_bases = 1;
        approximation.method = 'SSGP';
        nfeatures = 32;
        approximation.nfeatures = nfeatures;
        [approximation.phi, approximation.dphi_dx] = sample_features_GP(theta.cov, model, approximation);
        for i =1:m
            gs = sample_binary_GP_precomputed_features(xtrain, ctrain, theta.cov, model, approximation, post);
            fx_decoupled_Fourier(i, :)=gs(xtest);
        end
        
        W_decoupled_Fourier(k,j) = Wasserstein2(mu_f, Sigma2_f, fx_decoupled_Fourier);
    end
end

mr = 1;
mc = 3;
legend_pos = [-0.18,1];

fig=figure('units','centimeters','outerposition',1+[0 0 fwidth fheight(1)]);
fig.Color =  [1 1 1];
layout = tiledlayout(mr,mc, 'TileSpacing', 'tight', 'padding','tight');
i = 0;

nexttile();
i=i+1;
h1 = plot(num_tr, log10(W_decoupled_Hilbert), 'Color',  C(1,:),'LineWidth', linewidth); hold on;
h2 = plot(num_tr, log10(W_decoupled_Fourier), 'Color',  C(2,:),'LineWidth', linewidth); hold on;
h3= plot(num_tr, log10(W_weight_Hilbert), 'Color',  C(3,:),'LineWidth', linewidth); hold on;

box off
ylabel('log$_{10}$2-Wasserstein')
xlabel('Size of the training set')
% yl = get(gca,'Ylim');

set(gca, 'Xtick', num_tr,  'Fontsize', Fontsize);
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
legend([h1, h2, h3], 'Decoupled bases, RRGP', 'Decoupled bases, SSGP', 'Weight-space approximation, RRGP')
box off
legend box off
% figname  = 'GP_sampling_binary';
% folder = [figure_path,figname];
% savefig(fig, [folder,'/', figname, '.fig'])
% exportgraphics(fig, [folder,'/' , figname, '.pdf']);
% exportgraphics(fig, [folder,'/' , figname, '.png'], 'Resolution', 300);
%

