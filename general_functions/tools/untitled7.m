x = linspace(0,1,100);
y = g(x);

figure();
plot(x, y)

[d,n]= size(x);

kernelname = base_name;
[~,  g_mu_y, g_sigma2_y, g_Sigma2_y, dmuc_dx, dmuy_dx, dsigmay_dx] = prediction_bin_preference(theta, xtrain_norm, ctrain, [x;x0*ones(1,n)], kernelfun,kernelname, 'modeltype', modeltype);

figure();
plot(x, g_mu_y);


test = results; 

figure();
plot(test.score)
figure();
plot(results.score)