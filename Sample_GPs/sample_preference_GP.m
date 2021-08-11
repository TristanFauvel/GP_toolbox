function [sample_f, sample_g, decomposition]= sample_preference_GP(x, theta, xtrain, ctrain, model, approximation, post)

D = size(xtrain,1)/2; %dimension

[phi_pref, dphi_pref_dx, phi, dphi_dx] = sample_features_preference_GP(theta, D, model, approximation);

    
[sample_g, dsample_g_dx, decomposition] = sample_value_GP_precomputed_features(phi, dphi_dx, phi_pref, dphi_pref_dx, xtrain, ctrain, theta, model, approximation, post);

 
sample_g  = sample_g(x);
sample_f = sample_g - sample_g';
sample_f = sample_f(:);

return

