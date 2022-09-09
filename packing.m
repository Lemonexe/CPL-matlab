%PACKING
%       With no argument supplied and no output, all packings defined in the table 'data/packing.csv'
%       When a number is supplied, a vector with values corresponding to packing index is returned into output
%       Taken from the appendix of "Liquid-Liquid and Solid-Liquid Extractors" (2016) by Jean-Paul Duroudier
%
%CALL EXAMPLE:
%       packing
%       data = packing(12)

function varargout = packing(varargin)
    addpath('lib');

    data = readCSVasCell('data/packing.csv');

    % return a specific packing
    if nargin == 1 && nargout == 1
        dim = size(data);
        pointer = varargin{1};
        if pointer > dim(1)
            error('Packing %d not found, database only contains packings from 1 to %d', pointer, dim(1));
        end
        varargout{1} = data(pointer, :);
    % print a table of all packings
    else
        a = size(data);
        fprintf(1, '%6s%18s%17s%23s%27s%26s%14s\n', 'index', 'type', 'dimension [mm]', 'wall thickness [mm]', 'apparent density [kg/m3]', 'volumetric surface [m–1]', 'porosity [%]');
        for i = 1 : a(1)
            if ~isnumeric(data{i, 4}); data{i, 4} = NaN; end
            fprintf(1, '%6i%18s%17.1f%23.1f%27.0f%26.0f%14.0f\n', i, [data{i, 1}, ' ', data{i, 2}], data{i, 3}, data{i, 4}, data{i, 5}, data{i, 6}, data{i, 7})
        end
    end
end
