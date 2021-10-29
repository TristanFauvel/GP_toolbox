function [sample_f, sample_g, decomposition]= sample_preference_GP(x, theta, xtrain_norm, ctrain, model, approximation, post)

D = model.D; %dimension

[approximation.phi_pref, approximation.dphi_pref_dx, approximation.phi, approximation.dphi_dx] = sample_features_preference_GP(theta, D, model, approximation);

    
[sample_g, dsample_g_dx, decomposition] = sample_value_GP_precomputed_features(approximation, theta, xtrain_norm, ctrain, model, post);

 
sample_g  = sample_g(x);
sample_f = sample_g - sample_g';
sample_f = sample_f(:);

return

