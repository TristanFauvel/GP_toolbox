function save_benchmark_results(acquisition_name, xtrain, xtrain_norm, ctrain, score, xbest, g, objective, pathname)

fi = ['xtrain_',acquisition_name];
experiment.(fi) = xtrain;
fi = ['xtrain_norm_',acquisition_name];
experiment.(fi) = xtrain_norm;
fi = ['ctrain_',acquisition_name];
experiment.(fi) = ctrain;
fi = ['score_',acquisition_name];
experiment.(fi) = score;
fi = ['xbest_',acquisition_name];
experiment.(fi) = xbest;
if save_data
    filename = [pathname,'/Data_cPBO/',objective,'_',acquisition_name];
    close all
    save(filename, 'experiment')
end