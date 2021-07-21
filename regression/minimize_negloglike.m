function [negL, dnegL] =minimize_negloglike(hyp, x_tr, y_tr, kernelfun, meanfun,ncov_hyp, nmean_hyp, update)
hyp = hyp(:);
if strcmp(update, 'none')
    negL = [];
    dnegL = zeros(1,ncov_hyp+nmean_hyp);
else
    theta.cov = hyp(1:ncov_hyp);
    theta.mean = hyp((ncov_hyp+1):(ncov_hyp+nmean_hyp));
    [negL, dnegL] = negloglike(theta, x_tr, y_tr, kernelfun, meanfun);
    if strcmp(update, 'cov')
        dnegL(ncov_hyp+1:end)=0;
    elseif strcmp(update, 'mean')
        dnegL(1:ncov_hyp)=0;
    end
end

return


