function [phi_pref, dphi_pref_dx, phi, dphi_dx] = sample_features_preference_GP(theta, model, approximation)
%% SSGP : Method based on the Sparse-Spectrum GP, Lazaro-Gredilla 2010
%% RRGP: Method based on the Reduced-Rank GP, Solin 2019
[phi, dphi_dx] = sample_features_GP(theta, model, approximation);

phi_pref = @(x)  phi(x(1:model.D,:)) -  phi(x(model.D+1:end,:));

dphi_pref_dx = @(x)  dphi_dx(x(1:model.D,:));
return
