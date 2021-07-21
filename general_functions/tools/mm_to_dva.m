function rdva = mm_to_dva(r)
% Formula from Watson, Journal of Vision 2014
rdva = 3.556*r+0.05993*r.^2-0.007358*r.^3+0.0003027*r.^4;
