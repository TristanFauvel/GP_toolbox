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
     filename = [pathname,objective,'_',acquisition_name];
    save(filename, 'experiment')
 