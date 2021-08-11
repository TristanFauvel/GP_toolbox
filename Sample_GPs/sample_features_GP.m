function [phi, dphi_dx] = sample_features_GP(theta, D, model, approximation)
%% SSGP : Method based on the Sparse-Spectrum GP, Lazaro-Gredilla 2010
%% RRGP: Method based on the Reduced-Rank GP, Solin 2019

nfeatures = approximation.nfeatures;
kernelname = model.kernelname;
switch approximation.method
    case 'SSGP'

 %       m = 6561; %as m increases, the approximation gets better
        %only valid for Gaussian_kernelfun_wnoise
        % lambda       =  exp(theta(1));
        % k0           =  exp(theta(2));

        if strcmp(kernelname, 'Gaussian_wnoise')
            lambda       =  exp(theta(1))*ones(1,D);
            k0           =  exp(theta(2));
        elseif strcmp(kernelname, 'Gaussian')
            lambda       =  exp(theta(1))*ones(1,D);
            k0           =  exp(theta(2));
        elseif strcmp(kernelname, 'ARD')
            lambda = exp(theta(1:end-1));
            k0 = exp(theta(end));
        elseif strcmp(kernelname, 'ARD_wnoise')
            lambda = exp(theta(1:end-2));
            k0 = exp(theta(end-1));
%         elseif strcmp(kernelname, 'Matern52')
%             lambda = exp(theta(1));
%             k0 = exp(theta(2));
        else
            error('The SSGP approximation.method is not implemented for this kernel')
        end
        sigma = diag(lambda)/((2*pi)^2);

        mu = zeros(1,D);
        W= mvnrnd(mu, sigma, nfeatures); % W is sampled from the normalized spectral density of the gaussian kernel.
        b = 2 * pi * rand(nfeatures, 1); %b is sampled uniformly in [0, 2pi]

        alpha = k0;
        phi = @(x) sqrt(2*alpha/nfeatures)*cos(W*x*(2*pi)+b)'; % phi : ntest x nfeatures
        dphi_dx = @(x) (-2*pi*sqrt(2*alpha/nfeatures)*sin((2*pi)*W*x+b).*W); % nfeatures x D
    case 'RRGP' %Reduced-rank approximation.method, Solin 2019
        
        if numel(theta) ~=2
           error('The number of hyperparameters for Matérn kernels is 2') 
        end
        if D ~= 8
            m =floor(nfeatures^(1/D)); %/D;6561
            j= myndgrid(1:m,D);
        else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nfeatures =  [3, 3, 3, 3, 3, 3, 3, 3];%
            j = ndgrid_ndinputs(nfeatures);
        end
        nfeatures = size(j,1);
        L = 2; %1
        switch kernelname
            case 'Matern32'
                nu = 3/2; %or 5/2 depending of the Matérn kernel
            case 'Matern52'
                nu =5/2; %or 5/2 depending of the Matérn kernel
        end

        sqrtlambda = sqrt(sum((pi*j./(2*L)).^2,2)); %(61) In Solin 2019

        if D == 1
            phi_j = @(x) sin(sqrtlambda.*(x+L))./sqrt(L);%(60) In Solin 2019
            dphi_j_dx = @(x) (sqrtlambda.*cos(sqrtlambda.*(x+L))./sqrt(L));
        else
            jj = @(x) permute( repmat(j,1,1,size(x,2)), [2,3,1]);
            arg  = @(x) pi*(jj(x).*(x+L))/(2*L);
            phi_j = @(x) permute(prod(sin(arg(x))./sqrt(L),1), [3,2,1]);%(60) In Solin 2019

            dphi_j_dx = @(x) feature_derivative(x,L, jj, D, nfeatures, arg);

        end
        if strcmp(kernelname, 'Matern32') || strcmp(kernelname, 'Matern52')
        sqrtSlambda = sqrt(Matern_spectral_density(sqrtlambda, nu, theta, D));
        elseif strcmp(kernelname, 'ARD')
             sqrtSlambda = sqrt(SE_ARD_spectral_density(sqrtlambda, theta, D));
        end
        phi = @(x) (sqrtSlambda.* phi_j(x))'; % ntest x nfeatures
        %             dphi_dx = @(x) squeeze((sqrtSlambda.*ones(nfeatures,size(x,2))).* dphi_j_dx(x))';
        dphi_dx = @(x) (sqrtSlambda.*ones(nfeatures,size(x,2))).* squeeze(dphi_j_dx(x)); % nfeatures x D
end

return


function dphi_j_dx = feature_derivative(x,L, jj, D, nfeatures, arg)

sinvals = sin(arg(x))./sqrt(L);
cosvals = pi*jj(x).*cos(arg(x))./(2*L^1.5);

dphi_j_dx = zeros(size(x,2),nfeatures,D);
for d = 1:D
    dims = [1:d-1,d+1:D];
    dphi_j_dx(:,:,d) = squeeze(prod([sinvals(dims,:,:);cosvals(d,:,:)],1));
end
