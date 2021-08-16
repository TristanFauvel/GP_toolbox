function [sample_f, dsample_f_dx] = sample_GP(theta, xtrain, ytrain, model, approximation)

D = size(xtrain,1);
[phi,dphi_dx] = sample_features_GP(theta, model, approximation);
[sample_f, dsample_f_dx] = sample_GP_precomputed_features(theta, phi, dphi_dx, xtrain, ytrain(:), model, approximation);
