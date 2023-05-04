function varargout = beautifulbox(data, groups, pointopt, boxopt, varargin)
% h = BEAUTIFULBOX(data, categories, [], [], boxplotoption);
% h = BEAUTIFULBOX(data, categories, [], boxopt, boxplotoption);
% h = BEAUTIFULBOX(data, categories, pointopt, boxopt, boxplotoption);
% h = BEAUTIFULBOX(data, categories, pointopt, boxopt, boxplotoption);
%
% Plot a box plot with built-in boxplot function and plot corresponding
% datapoints.
%
% INPUTS :
% ------
% data : data to box plot.
% categories : array that describes categories of each data point.
% pointopt : struct with parameters for a line object, eg. options.LineWidth
% = 0.1, options.Marker = 'o', plus a field 'jitter' giving a random x axis
% offset for data points. Set MarkerFaceColor property to 'same' to get the
% same color than the corresponding box.
% boxopt : struct, tune rendering of box plot. See boxchart doc to see
% available options. A field BoxColors can be set a color for each box, in
% this case it has to be ngroupsx3 array. The field Labels will be used to
% label each groups.
% Everything else is sent to the boxplot function after data and groups.

% RETURNS :
% -------
% output of the boxplot function.

% --- Check input
p = inputParser;
p.addRequired('data', @isnumeric);
p.addRequired('groups', @(x) isnumeric(x)||ischar(x)||isstring(x)||iscell(x)||iscategorical(x));
p.addRequired('pointopt', @(x) isempty(x)||isstruct(x));
p.addRequired('boxopt', @(x) isempty(x)||isstruct(x));
p.parse(data, groups, pointopt, boxopt);

data = p.Results.data;
groups = p.Results.groups;
popt = p.Results.pointopt;
bopt = p.Results.boxopt;

g = unique(groups);

if isempty(popt)
    popt = struct;
end
if isempty(bopt)
    bopt = struct;
end
if isfield(popt, 'jitter')
    jitter = popt.jitter;
    popt = rmfield(popt, 'jitter');
else
    jitter = 0;
end
if isfield(popt, 'Color')
    pointscolors = popt.Color;
    popt = rmfield(popt, 'Color');
else
    pointscolors = lines(numel(g));
end
if isfield(bopt, 'BoxColors')
    boxcolors = bopt.BoxColors;
    bopt = rmfield(bopt, 'BoxColors');
end
if isfield(bopt, 'Labels')
    boxlabels = bopt.Labels;
    bopt = rmfield(bopt, 'Labels');
end

% --- Defaults values

% - Data points
defpopt = struct;
defpopt.LineWidth = 1;
defpopt.LineStyle = 'none';
defpopt.Marker = 'o';
defpopt.MarkerFaceColor = 'none';
defpopt.MarkerSize = 6;
opt_points = addToStructure(defpopt, popt);

% - Boxes
defbopt = struct;
defbopt.LineWidth = 1.25;
defbopt.MarkerStyle = '+';
defbopt.MarkerColor = 'r';
defbopt.MarkerSize = 8;
opt_boxes = addToStructure(defbopt, bopt);

% --- Create figure
ax = gca; hold on

% --- Add data points & box
for idx = 1:numel(g)
    
    datapoints = data(ismember(groups, g(idx)));
    x = idx.*ones(1, numel(datapoints));
    xjitter = x + rescale(rand([1, numel(datapoints)]), -jitter/2, jitter/2);
    p = plot(ax, x, datapoints);
    p.Color = pointscolors(idx, :);
    
    % Check options
    if isfield(pointopt, 'MarkerFaceColor') && strcmp(pointopt.MarkerFaceColor, 'same')
        opt_points.MarkerFaceColor = p.Color;
    end
    
    p = mergeStructures(opt_points, p);
    
    if exist('boxcolors', 'var')
        b = boxchart(ax, x, datapoints, 'BoxFaceColor', boxcolors(idx, :), varargin{:});
    else
        b = boxchart(ax, x, datapoints, varargin{:});
    end
    
    b = mergeStructures(opt_boxes, b);
    
    % Apply jitter
    p.XData = xjitter;
end

% Style adjustment
ax.XTick = 1:numel(g);
ax.XTickLabel = boxlabels;
ax.Box = 'on';

end