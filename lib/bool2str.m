%s = BOOL2STR(b)
%       converts logical value b into a string (to be used for display)

function s = bool2str(b)
    strs = {'false', 'true'};
    s = strs{b+1};
end
