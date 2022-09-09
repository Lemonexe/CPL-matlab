%data = READCSVASCELL(fileName)
%       reads CSV file to a 2D cell of strings or numbers (they are autodetected)
%       better than csvread, which reads the file as a numeric matrix, and cannot handle strings
%       supply 'fileName' of the .csv file

function data = readCSVasCell(fileName)
    fileID = fopen(fileName, 'r');

    % init
    data = {};
    i = 1;
    tline = fgetl(fileID);

    % read whole file with a while cycle one line at a time
    while tline >= 0
        row = split(tline, ',');

        % foreach cell: if number, assign as number, else as a string
        for j = 1:numel(row)
            [X, tf] = str2num(row{j});
            if tf; data{i, j} = X;
            else data{i, j} = row{j};
            end
        end
        i = i + 1;
        tline = fgetl(fileID);
    end
    fclose(fileID);
    %that all could have been expressed by a single line, which would, however, require Matlab 2013:
    %data = table2cell(readtable(fileName))
end
