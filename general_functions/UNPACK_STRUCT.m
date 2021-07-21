function UNPACK_STRUCT(strct, warn)
%UNPACK_STRUCT make all fields of the structure separate variables in the calling workspace
%
%     UNPACK_STRUCT(strct[, warn])
%
% All of the fields of strct will now exist as top-level variables in the
% calling workspace. This may remove a lot of 'strct.' clutter from code.
%
% Unless the warn option is given and set to 'false', the code will emit
% warnings if pre-existing variables are being over-written.

% Iain Murray, October 2009

%DEFAULT('warn', true);

args = fieldnames(strct);

for ff = args(:)'
    field = ff{1};
    if warn && evalin('caller', ['exist(''', field, ''', ''var'')'])
        warning(['Over-writing variable ''', field, '''. Set warn=false if this was intended']);
    end
    assignin('caller', field, strct.(field));
end