%[w, u] = GETEQFUN(sol, SA, t, w_aq)
%       Obtains tabelated equilibrium function for a given temperature,
%
%PARAMETERS:
%       sol = organic solvent = 'TOL', 'TCE'
%       SA = content of ammonium sulphate % = 30, 35, 40 etc.
%       t_arg = requested temperature in °C, warning will be thrown if limits are exceeded
%       OPTIONAL aq = vector of aqueous mass frac values where the function will be tabulated
%
%OUTPUT varies depending on number of requested output arguments: (w is mass frac in aquous phase, u in organic)
%       - if there are 0, plot is drawn as u(w)
%       - if there are 2, tabelated function as [w, u] is returned
%       - if there are 4, tabelated function as [w, u, max(w), max(u)] is returned
%
%CALL EXAMPLE:
%       getEqFun('TOL', 40, 60)
%       [w, u] = getEqFun('TOL', 40, 60)

function varargout = getEqFun(varargin)
    addpath('lib'); % dependencies

    if nargin >= 3; sol = varargin{1}; SA = varargin{2}; t_arg = varargin{3}; end
    if nargin == 4; aq = varargin{4}; else; aq = 0; end

    % load all datasets from file and convert the result into a flat 2D cell
    data = readCSVasCell('data/params.csv');

    % find parameters
    rowName = sprintf('%s_SA%d', sol, SA);
    found = false;
    m = size(data); m = m(1);
    for i = 1:m
        if strcmp(rowName, data{i, 1});
            found = true;
            t_min  = data{i, 2};
            t_max  = data{i, 3};
            aq_max = data{i, 4};
            params = cell2mat(data(i, 5:end));
            break;
        end
    end
    if ~found; error('Data for %s SA %d %% not found!', sol, SA); end

    % warn if required temperature is out of bounds. Do not track warning for brevity & clarity, as it's intended for end user and not for development
    if t_arg < t_min || t_arg > t_max
        warning('off', 'backtrace'); warning('EXTRAPOLATION\ntemperature (%d°C) should be within bounds %d-%d°C', t_arg, t_min, t_max); warning('on', 'backtrace');
    end

    % equilibrium function tabelation
    aux = size(aq);
    % if undefined aq vector, elseif defined, but needs to be transposed to column vector
    if length(aq) == 1
        step = 1e-5;
        aq = (0 : step : (aq_max + step))';
    elseif aux(1) == 1
        aq = aq';
    end
    t = ones(size(aq)) * t_arg;
    org = eqFunT(params, [t aq]);

    % two alternatives for function OUTPUT
    % IF there are two output parameters defined, then fill them with tabelation values
    if nargout == 2
        varargout = {aq, org};
    elseif nargout == 4
        varargout = {aq, org, aq_max, eqFunT(params, [t_arg aq_max])};
    else
        % ELSE draw plots of original data and the tabelated function
        plotLLEdata({'big','title'}, sol, SA);
        hold on;
        plot(aq*100, org*100, '-');
        hold off;
    end
end
