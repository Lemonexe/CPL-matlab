%nbytes = HEADLINE(s)
%       superstructure for fprintf - prints input string and an ascii art underline made of diaereses

function nbytes = headline(s)
    nbytes = fprintf(1, [s, '\n', repmat('¨', 1, length(s)), '\n']);
end
