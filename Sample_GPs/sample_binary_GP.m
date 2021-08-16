function [sample_g, dsample_g_dx, decomposition] = sample_binary_GP(theta, xtrain, ctrain, model, approximation, post)
%% SSGP : Method based on the Sparse-Spectrum GP, Lazaro-Gredilla 2010
%% RRGP: Method based on the Reduced-Rank GP, Solin 2019
% x_data : (2*dimension)*n
%y_data : n x 1

[phi, dphi_dx] = sample_features_GP(theta, model, approximation);
[sample_g, dsample_g_dx, decomposition] = sample_binary_GP_precomputed_features(phi, dphi_dx, xtrain, ctrain, theta,model, approximation, post);

return