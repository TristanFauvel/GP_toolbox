function save_benchmark_results(acquisition_name, structure_name, xtrain, ctrain, score, xbest, objective, pathname, task)

fi = ['xtrain_',structure_name];
experiment.(fi) = xtrain;
fi = ['ctrain_',structure_name];
experiment.(fi) = ctrain;
fi = ['score_',structure_name];
experiment.(fi) = score;
fi = ['xbest_',structure_name];
experiment.(fi) = xbest;
filename = [pathname,'/',task, '_', objective,'_',acquisition_name,'.mat'];

if ~isfolder(pathname)
    mkdir(pathname)
end
save(filename, 'experiment')
