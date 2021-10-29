function [sample_f, dsample_f_dx] = sample_GP(theta, xtrain, ytrain, model, approximation)

[approximation.phi, approximation.dphi_dx] = sample_features_GP(theta, model, approximation);
[sample_f, dsample_f_dx] = sample_GP_precomputed_features(xtrain, ytrain(:), theta, model, approximation);

