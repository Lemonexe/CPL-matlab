%arr = SPLIT(str, delim)
%       splits string into a cell array of substrings delimited by delimiter

function arr = split(str, delim)
    dels = strfind(str, delim); % positions of all delimiters
    if isempty(dels); arr = {str}; return; end % no delimiters, return single cell with entire string
    arr = cell(1, length(dels)+1);
    arr{1,1} = str(1:(dels(1)-1)); % first string
    for i = 2:length(dels) % all between delimiters
        arr{1,i} = str((dels(i-1)+1):(dels(i)-1));
    end
    arr{1,end} = str((dels(end)+1):end); % last string
end
