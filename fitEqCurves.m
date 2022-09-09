%FITEQCURVES
%       one-time routine which reads all CSV raw data files in the folder data/LLE/ and performs curve fitting with equilibrium function 'eqFunT'
%       results (parameters from regression) are saved in data/params.csv

function fitEqCurves
    addpath('lib'); % dependencies

    % read contents of data/ directory
    list = dir('data/LLE');
    list = {list.name};
    list = list(~ismember(list, {'.','..'}));
    M = [];

    % options for lsqcurvefit - setting lower TolFun is necessary to get a good fit
    options = optimset('TolFun', 1e-9, 'Display', 'off', 'MaxFunEvals', 3000, 'MaxIter', 1000);

    % go through each file in data/: read the file, call lsqcurvefit on it and save both data scope and fitted parameters
    for i = 1:length(list)
        file = csvread(['data/LLE/' list{i}]);
        xdata = file(:, 1:2);
        ydata = file(:, 3);

        x_max = max(xdata);
        x_min = min(xdata);
        params = lsqcurvefit(@eqFunT, ones(1,6), xdata, ydata, [], [], options);

        % note that the first three values are: min temperature, max temperature, max CPL frac in aqueous phase (min CPL is always zero)
        M(i, :) = [x_min(1), x_max, params];
    end

    % generate content of the new CSV row by row
    fileID = fopen('data/params.csv', 'w');
    for i = 1:length(list)
        fprintf(fileID, [strrep(list{i}, '.csv', '') ',' joinVector(M(i, :), ',') '\n']);
    end
    fclose(fileID);

    fprintf(1, 'Success!\n');
end

% like strjoin, but joins a numeric vector
function str = joinVector(vec, delim)
    str = '';
    for i = 1:length(vec)
        str = [str sprintf('%.10g', vec(i))];
        if i < length(vec); str = [str delim]; end
    end
end
