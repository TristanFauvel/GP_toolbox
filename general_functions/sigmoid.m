function y= sigmoid(x, x0, k)
DEFAULT('x0',0)
DEFAULT('k', 1)
y = 1 ./ (1 + exp(-k*(x-x0)));
return