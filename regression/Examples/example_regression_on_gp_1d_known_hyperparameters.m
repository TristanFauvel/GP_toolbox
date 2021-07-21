clear all;
close all;

addpath(genpath('../GP_toolbox'))

n=1000;
rng(10)

meanfun= @constant_mean;

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

theta_true = [1,2];


x = linspace(0,5,n);

y = mvnrnd(meanfun(x,0), kernelfun(theta_true, x,x));

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
    x_tr = x(:,i_tr);
    y_tr = y(:, i_tr); %+ 0.1*randn(1,ntr);
    
    x_test = x;
    y_test = y;
    
    
    %     theta.cov = rand(ncov_hyp,1);
    %     theta.mean = zeros(nmean_hyp,1);
    theta.cov = theta_true;
    theta.mean = 0;
    
    [mu_y, sigma2_y]= prediction(theta, x_tr, y_tr, x_test, kernelfun, meanfun);
    
    h=figure(2);
    h.Name = 'Bayesian optimisation before hyperparameters optimization';
    h.Color =  [1 1 1];
    plot(x,y,'LineWidth',1.5); hold on;
    scatter(x_tr, y_tr) ; hold on;
    errorshaded(x, mu_y, sqrt(sigma2_y), 'Color', 'red','DisplayName','Prediction'); hold off
    xlabel('x','Fontsize',14)
    ylabel('f(x)','Fontsize',14)
    grid off
    
        
    %% Prediction with the new hyperparameters
    [mu_y, sigma2_y]= prediction(theta, x_tr, y_tr, x_test, kernelfun, meanfun);
    
    
    c= othercolor('GnBu7');
    K= kernelfun(theta_true, x(1:10:end),x(1:10:end));
    
    xplot =  x(1:10:end);
    [~,~,~, ~, K] =  prediction(theta, x_tr, y_tr, xplot, kernelfun, meanfun);
    
    h=figure(2);
    h.Name = 'Bayesian optimisation after hyperparameters optimization';
    h.Color =  [1 1 1];
    plot(x,y,'LineWidth',1.5); hold on;
    scatter(x_tr, y_tr) ; hold on;
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
    scatter(x_tr, y_tr) ; hold on;
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