function argStruct = namevaluepairtostruct(defaults, varargin)
%NAMEVALUEPAIRTOSTRUCT Converts name/value pairs to a struct.
% 
% ARGSTRUCT = NAMEVALUEPAIRTOSTRUCT(DEFAULTS, VARARGIN) converts
% name/value pairs to a struct, with defaults.  The function expects an
% even number of arguments to VARARGIN, alternating NAME then VALUE.
% (Each NAME should be a valid variable name.)
% 
% Examples: 
% 
% No defaults
% NameValuePairToStruct(struct, ...
%    'foo', 123, ...
%    'bar', 'qwerty', ...
%    'baz', magic(3))
% 
% With defaults
% NameValuePairToStruct( ...
%    struct('bar', 'dvorak', 'quux', eye(3)), ...
%    'foo', 123, ...
%    'bar', 'qwerty', ...
%    'baz', magic(3))
% 
% See also: inputParser
varargin = varargin{:};
nArgs = length(varargin);
if rem(nArgs, 2) ~= 0
   error('NameValuePairToStruct:NotNameValuePairs', ...
      'Inputs were not name/value pairs');
end

argStruct = defaults;
for i = 1:2:nArgs
   name = varargin{i};
   if ~isvarname(name)
      error('NameValuePairToStruct:InvalidName', ...
         'A variable name was not valid');
   end
   argStruct = setfield(argStruct, name, varargin{i + 1});  %#ok<SFLD>
end

end