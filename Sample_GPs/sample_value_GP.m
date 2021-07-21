function [sample_g, dsample_g_dx, decomposition] = sample_value_GP(theta, xtrain, ctrain, base_kernelname, kernelfun, decoupled_bases, modeltype, approximation_method, base_kernelfun, nfeatures, post, condition)
%% SSGP : Method based on the Sparse-Spectrum GP, Lazaro-Gredilla 2010
%% RRGP: Method based on the Reduced-Rank GP, Solin 2019
D = size(xtrain,1)/2; %dimension

% xtrain : (2*dimension)*n
%y_data : n x 1

[phi_pref, dphi_pref_dx, phi, dphi_dx] = sample_features_preference_GP(theta, D, base_kernelname, approximation_method, nfeatures);
[sample_g, dsample_g_dx, decomposition] = sample_value_GP_precomputed_features(phi, dphi_dx, phi_pref, dphi_pref_dx, xtrain, ctrain, theta,kernelfun, decoupled_bases, modeltype, base_kernelfun, post, condition);

return

