%PLOTLLEDATA
%       this function will generate plots of measured data using the CSV files in the folder data/LLE/
%
%PARAMETERS are all optional:
%       opts = output options as cell of strings, can contain following strings to trigger actions:
%              KD    = also draw plots of distribution coefficient
%              save  = also save plots as .png image files
%              big   = make plots fullscreen
%              gray  = make black and white by using different markers instead of using colors are used
%              title = write title of plot
%       sol  = only the selected dataset will be displayed (undefined by default)
%       SA   = used together with sol (undefined by default)
%
%OUTPUT
%       if 'save' parameter, plot images will be saved in output/ folder
%       if there are two arguments supplied as (solvent, SA), only the selected dataset will be plotted
%
%CALL EXAMPLE:
%       plotLLEdata(false,  'TOL', 40);

function plotLLEdata(varargin)
    optsCell = false;
    sol_arg = false; % false or string
    SA_arg = false; % false or string or number
    found = false; % found match for sol and SA
    if nargin >= 1; optsCell = varargin{1}; end
    if nargin == 3; sol_arg = varargin{2}; SA_arg = varargin{3}; end
    if nargin > 3; error('function plotLLEdata can only accept up to three parameters, but %d encountered!', nargin); end

    if ~iscell(optsCell); optsCell = {}; end
    opts = struct();
    for i = 1:numel(optsCell); opts.(optsCell{i}) = true; end

    opts.KD   = isfield(opts, 'KD')   && opts.KD;
    opts.save = isfield(opts, 'save') && opts.save;
    opts.big  = isfield(opts, 'big')  && opts.big;
    opts.gray = isfield(opts, 'gray') && opts.gray;
    opts.title= isfield(opts, 'title')&& opts.title;
    if opts.KD && opts.gray; warning('distibution coefficient plots cannot be drawn in grayscale'); end

    if opts.save && ~exist('output', 'dir'); mkdir('output'); end

    list = dir('data/LLE');
    list = {list.name};
    list = list(~ismember(list, {'.','..'}));

    % draw a plot: XY means whether to draw u(w), else KD(w)
    function doPlot(XY)
        figure();
        % different colors to be used
        colors = {'r', 'm', 'g', 'b', 'k', 'y', 'c'};
        % or different markers to be used
        markers = {'o', 's', '^', 'v', 'd', 'h', 'p', 'x', '+'};

        hold on;
        for j = 1:numel(datasets)
            s = datasets{j} * 100; % to percentage
            s = sortrows(s, 1);

            if opts.gray && XY
                plot(s(:, 1), s(:, 2), 'Marker', markers{j}, 'MarkerFace', 'none', 'MarkerEdge', 'k', 'LineStyle', 'none');
            elseif XY
                plot(s(:, 1), s(:, 2), ['x' colors{j}]);
            else
                plot(s(:, 1), s(:, 2) ./ s(:, 1), ['-' colors{j}]);
            end
        end
        hold off;
        legobj = legend(legends, 'Location', 'eastoutside');
        legobj.Box = 'off';
        xlabel('\itw\rm [%]');
        if XY; ylabel('\itu\rm [%]'); else; ylabel('\itK_{D}\rm'); end
        if opts.title; title(sprintf('%s, SA %s %%', sol, SA)); end
        if opts.big; set(gcf, 'Position', get(0, 'Screensize')); else; set(gcf, 'Position', [100, 100, 600, 400]); end
        box on;

        % save figure 'filename' in 'dirname'
        if opts.save;
            if XY; dirname = 'output/XY'; else dirname = 'output/KD'; end
            filename = sprintf('%s_SA%s.png', sol, SA);
            if ~exist(dirname, 'dir'); mkdir(dirname); end
            print([dirname '/' filename],'-dpng', '-r300');
         end
    end

    % all files in data/
    for i = 1:numel(list)
        % get solvent and SA from filename
        [s, e] = regexp(list{i}, '^[^_]+');
        sol = list{i}(s:e);

        [s, e] = regexp(list{i}, 'SA\d+');
        SA = list{i}((s+2):e);

        if nargin == 3 && (~strcmp(sol, sol_arg) || ~strcmp(SA, num2str(SA_arg))); continue; end
        found = true;

        % read data and plot it
        file = csvread(['data/LLE/' list{i}]);
        t   = file(:, 1); % [°C]
        aq  = file(:, 2); % absolute frac
        org = file(:, 3); % absolute frac

        t_dis = sort(unique(t), 'descend'); % distinct temperature values

        datasets = cell(1, length(t_dis));
        legends = cell(1, length(t_dis));
        for j = 1:length(t_dis)
            datasets{j} = [aq(t == t_dis(j)), org(t == t_dis(j))];
            legends{j} = sprintf('%d°C', t_dis(j));
        end

        doPlot(true);
        if opts.KD; doPlot(false); end
    end

    if nargin == 3 && ~found; error('File for %s, SA %d %% not found!', sol_arg, SA_arg); end
end
