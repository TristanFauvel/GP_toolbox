clear all;
close all;

addpath(genpath('../GP_toolbox'))

n=1000;
rng(10)

meanfun= @constant_mean;
lb = 0; ub = 5;
% %ARD kernel with noise
% ncov_hyp=3;
% nmean_hyp=1;
% kernelfun = @ARD_kernelfun_wnoise;
% theta_true = [1,2,0];
%
%ARD kernel without noise
ncov_hyp=2;
nmean_hyp=1;
kernelfun = @ARD_kernelfun;
kernelfun = @Matern52_kernelfun;
kernelname = 'Matern52';
theta_true = [1,2];

D = 1;
regularization = 'nugget';
type = 'regression';
hyps.ncov_hyp =3; % number of hyperparameters for the covariance function
hyps.nmean_hyp =1; % number of hyperparameters for the mean function
hyps.hyp_lb = -10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
hyps.hyp_ub = 10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
model = gp_regression_model(D, meanfun, kernelfun, regularization, hyps, lb, ub, kernelname);


x = linspace(lb, ub, n);

y = mvnrnd(meanfun(x,0), kernelfun(theta_true, x,x, true, 'nugget'));

h=figure(1);
h.Name = 'True function';
h.Color =  [1 1 1];
plot(x, y,'LineWidth',1.5)
xlabel('x','Fontsize',14)
ylabel('f(x)','Fontsize',14)
grid off

maxiter = 20;
F(maxiter) = struct('cdata',[],'colormap',[]);
idx= randsample(n,maxiter);

for ntr= 1:maxiter
    i_tr = idx(1:ntr);
    xtrain = x(:,i_tr);
    y_tr = y(:, i_tr); %+ 0.1*randn(1,ntr);
    
    x_test = x;
    y_test = y;
    
    
    %     theta.cov = rand(ncov_hyp,1);
    %     theta.mean = zeros(nmean_hyp,1);
    theta.cov = theta_true;
    theta.mean = 0;
    
    [mu_y, sigma2_y]= model.prediction(theta, xtrain, y_tr, x_test,[]);
    
    h=figure(2);
    h.Name = 'Bayesian optimisation before hyperparameters optimization';
    h.Color =  [1 1 1];
    plot(x,y,'LineWidth',1.5); hold on;
    scatter(xtrain, y_tr) ; hold on;
    errorshaded(x, mu_y, sqrt(sigma2_y), 'Color', 'red','DisplayName','Prediction'); hold off
    xlabel('x','Fontsize',14)
    ylabel('f(x)','Fontsize',14)
    grid off
    
        
    %% Prediction with the new hyperparameters
    [mu_y, sigma2_y]= model.prediction(theta, xtrain, y_tr, x_test, []);
    
    
    c= othercolor('GnBu7');
    K= kernelfun(theta_true, x(1:10:end),x(1:10:end), [], regularization);
    
    xplot =  x(1:10:end);
    [~,~,~, ~, K] =  model.prediction(theta, xtrain, y_tr, xplot, []);
    
    h=figure(2);
    h.Name = 'Bayesian optimisation after hyperparameters optimization';
    h.Color =  [1 1 1];
    plot(x,y,'LineWidth',1.5); hold on;
    scatter(xtrain, y_tr) ; hold on;
    errorshaded(x, mu_y, sqrt(sigma2_y), 'Color', 'red','DisplayName','Prediction'); hold off
    xlabel('x','Fontsize',14)
    ylabel('f(x)','Fontsize',14)
    ylim([-5 5])
    box off 

    
       
    h=figure(3);
    set(h, 'Position', get(0, 'Screensize'));
    h.Name = 'Bayesian optimisation after hyperparameters optimization';
    h.Color =  [1 1 1];
    
    f2 =subplot(1,2,2);    
    imagesc(K)
    colormap(c)
    grid off
    set(gca,'visible','off')
    box off
    pbaspect([1 1 1])
    
    f1=subplot(1,2,1);
    h.Name = 'Bayesian optimisation after hyperparameters optimization';
    h.Color =  [1 1 1];
    plot(x,y,'LineWidth',1.5); hold on;
    scatter(xtrain, y_tr) ; hold on;
    errorshaded(x, mu_y, sqrt(sigma2_y), 'Color', 'red','DisplayName','Prediction'); hold off
    ylim([-5 5])
    box off 
    set(gca,'XColor', 'none','YColor','none')
    
    F(ntr) = getframe(gcf);

end

fps=2;
writerObj = VideoWriter('GP_regression.avi');
writerObj.FrameRate = fps;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);