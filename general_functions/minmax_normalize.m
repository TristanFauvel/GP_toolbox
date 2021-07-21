function [x_norm, min_x, max_x] = minmax_normalize(x) 
%% Normalization of the training set
min_x= min(x, [], 2);
max_x= max(x, [], 2);

x_norm = (x-min_x)./(max_x - min_x);
