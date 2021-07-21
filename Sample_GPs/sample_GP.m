function [sample_f, dsample_f_dx] = sample_GP(theta, xtrain, ytrain, kernelname,approximation_method, decoupled_bases, nfeatures, kernelfun)

D = size(xtrain,1);
[phi,dphi_dx] = sample_features_GP(theta, D, kernelname, approximation_method, nfeatures);
[sample_f, dsample_f_dx] = sample_GP_precomputed_features(theta, phi, dphi_dx, xtrain, ytrain(:), kernelname, decoupled_bases, kernelfun);

