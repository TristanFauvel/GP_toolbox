function varargout = mydeal(fun, inputs)

outputs = cell(1, nargout);
[outputs{:}] = fun(inputs{:});

varargout = outputs;
return
