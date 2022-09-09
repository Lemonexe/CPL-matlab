%[x, TP] = EXTRACTORFEED(v, sol, SA, t, varName, varRange)
%       for explanation of 'v', 'sol', 'SA' and 't', see function extractor.m
%       this function enables you to feed multiple values into extractor()
%       the modified variable is specified in 'varName' parameter, which can have following values: 'mF', 'mS', 'wF', 'wR', 'uS', 'uE', 'SA', 't'
%       all values of 'varRange' vector are fed into extractor() as the selected variable
%       if there are two output variables specified, the first is just varRange, the second is corresponding number of theoretical plates
%
%CALL EXAMPLE:
%       [x, TP] = extractorFeed([27, NaN, 0.01, 1e-4, 0, NaN], 'TOL', 40, 60, 'mS', [3.5 4:1:15]);

function varargout = extractorFeed(v, sol, SA, t, varName, varRange)
    % find index of the modified value in 'v' using the labels in varNames
    varNames = {'mF', 'mS', 'wF', 'wR', 'uS', 'uE', 'SA', 't'};
    index = find(strcmpi(varNames, varName));
    if isempty(index)
        fprintf(1, 'Input varName %s invalid. These are possible: %s \n', varName, strjoin(varNames, ', ')); return
    end

    % table headings
    fprintf(1, '\n\t%s, SA %d, %d°C\n', sol, SA, t);
    fprintf(1, '%s\tTP\tminS\n', varName);

    % for every value in varRange, perform extractor() simulation, save TP value and print the values
    TP = zeros(size(varRange));
    warning('off');
    for i = 1:length(varRange)
        if index <= 6; v(index) = varRange(i); % value is one of the balance variables
        elseif index == 7; SA = varRange(i); % value is content of SA - useless option really, as only 30,35,40 is allowed, but why not...
        elseif index == 8; t = varRange(i); % value is temperature
        end
        s = extractor(v, sol, SA, t);
        TP(i) = s.TP;
        fprintf(1, '%.2e\t%.2f\t%.2e\n', varRange(i), s.TP, s.mSmin);
    end
    warning('on');
    fprintf(1, '\n');

    % if output is required, feed the table of values & TP into it
    if nargout == 2
        varargout{1} = varRange;
        varargout{2} = TP;
    end
end
