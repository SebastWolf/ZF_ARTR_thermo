function [h, p] = errorshaded(x, y, e, varargin)
% h = ERRORSHADED(x, y, errors)
% h = ERRORSHADED(_, 'line', lineoptions)
% h = ERRORSHADED(_, 'patch', patchoptions)
% h = ERRORSHADED(_, 'ax', axobject)
%
% ERRORSHADED plot y versus x with errors e represented as a shaded area.
%
% INPUTS :
% ------
% x : 1D vector, x-axis.
% y : 1D vector, y-axis.
% e : 1D vector, errors, same size as x and y. If scalar, assumes errors 
% is the same for all points.
% 'line', options : structure, options for the line object, + 'LineAlpha'
% 'patch', options : structure, options for the line object.
% 'ax', options : axes object, defines axes where to plot
%
% RETURNS :
% -------
% h : handle to the line
% p : handle to the shaded area (patch object)

% --- Check input
p = inputParser;
p.addRequired('x', @isnumeric);
p.addRequired('y', @isnumeric);
p.addRequired('e', @isnumeric);
p.addParameter('line', [], @(x) isempty(x)||isstruct(x));
p.addParameter('patch', [], @(x) isempty(x)||isstruct(x));
p.addParameter('ax', gca, @(x) isa(x, 'matlab.graphics.axis.Axes'));
p.parse(x, y, e, varargin{:});

x = p.Results.x;
y = p.Results.y;
e = p.Results.e;
lineoptions = p.Results.line;
patchoptions = p.Results.patch;
ax = p.Results.ax;

% Check dimension
if size(x, 1) ~= 1
    x = x';
end
if size(y, 1) ~= 1
    y = y';
end
if size(e, 1) ~= 1
    e = e';
end
if numel(e) == 1
    err = e;
    e = err*ones(size(y));
end

% Fill defaults
if isempty(lineoptions)
    lineoptions = struct;
end
if isempty(patchoptions)
    patchoptions = struct;
end
if isfield(lineoptions, 'LineAlpha')
    linealpha = lineoptions.LineAlpha;
    lineoptions = rmfield(lineoptions, 'LineAlpha');
else
    linealpha = [];
end

% Defaults for line
defline = struct;
defline.LineWidth = 1;
defline.LineStyle = '-';
defline.Marker = '.';
defline.MarkerFaceColor = 'none';
defline.MarkerEdgeColor = 'auto';
defline.MarkerSize = 15;
optline = addToStructure(defline, lineoptions);

% Defaults for shaded area
defpatch = struct;
defpatch.FaceAlpha = 0.2;
defpatch.EdgeColor = 'none';
optpatch = addToStructure(defpatch, patchoptions);

% --- Create data to be plotted
ybot = y - e;
ytop = y + e;

% --- Create figure
hold(ax, 'on');

% Set colors
if ~isfield(optline, 'Color')
    defcolor = ax.ColorOrder(ax.ColorOrderIndex, :);
    optline.Color = defcolor;
end
patchcol = optline.Color;
optline.Color = [optline.Color, linealpha];

% --- Plot patch & line
p = patch(ax, [x fliplr(x)], [ytop fliplr(ybot)], patchcol);
p = mergeStructures(optpatch, p);

h = plot(ax, x, y);
h = mergeStructures(optline, h);

end