function [pdf, centers, edges] = computePDF(xi, values, varargin)
% Computes the empirical probability density function underlying the values at
% the points in edges. Edges are the edges of each bin, ie, the first it the
% left side of the first bin and the last the right side of the last bin. The
% vector corresponding to bins centers is returned along the pdf.
%
% INPUTS :
% ------
% xi : bins edges or centers at which the pdf will be computed. Set to
% empty to use automatic binning.
% values : data.
% 'method', value : method to estimade pdf : 'hist' (default, simple
% normalized histogram) or 'kde' (kernel smoothing density).
% 'mode', value : specify if xi are bin edges or bin centers. Default is
% 'centers'.
% 'param', cell : additional arguments for histcount or ksdensity
% functions.
%
% OUTPUTS :
% -------
% pdf : empirical probability density function.
% centers : bin centers, pdf should be plotted against this.

% --- Check input
p = inputParser;
p.addRequired('xi', @isnumeric);
p.addRequired('values', @isnumeric);
p.addParameter('method', 'hist', @(x) ischar(x)||isstring(x));
p.addParameter('mode', 'centers', @(x) ischar(x)||isstring(x));
p.addParameter('param', {}, @iscell);
p.parse(xi, values, varargin{:});

xi = p.Results.xi;
values = p.Results.values;
method = p.Results.method;
mode = p.Results.mode;
param = p.Results.param;

% Make sure that returns are in the same orientation as inputs.
if size(xi, 1) ~= 1
    flipbinvec = true;
    xi = xi';
else
    flipbinvec = false;
end
if size(values, 1) ~= 1
    flipvalues = true;
    values = values';
else
    flipvalues = false;
end

% --- Processing
switch method
    case 'hist'
        
        if isempty(xi)
            [pdf, edges] = histcounts(values, 'Normalization', 'pdf', param{:});
            centers = edges(1:end-1) + diff(edges) / 2;
        else
            switch mode
                case 'centers'
                    centers = xi;
                    d = diff(centers)/2;
                    edges = [centers(1) - d(1), centers(1:end-1) + d, centers(end) + d(end)];
                    edges(2:end) = edges(2:end)+eps(edges(2:end));
                case 'edges'
                    edges = xi;
                    centers = edges(1:end-1) + diff(edges)/2;
                otherwise
                    error('Choose ''centers'' or ''edges'' for bins.');
            end
            pdf = histcounts(values, edges, 'Normalization', 'pdf', param{:});
        end
        
    case 'kde'
        
        if isempty(xi)
            [pdf, centers] = ksdensity(values, param{:});
            d = diff(centers)/2;
            edges = [centers(1) - d(1), centers(1:end-1) + d, centers(end) + d(end)];
            edges(2:end) = edges(2:end)+eps(edges(2:end));
        else
            switch mode
                case 'centers'
                    centers = xi;
                    d = diff(centers)/2;
                    edges = [centers(1) - d(1), centers(1:end-1) + d, centers(end) + d(end)];
                    edges(2:end) = edges(2:end)+eps(edges(2:end));
                case 'edges'
                    edges = xi;
                    centers = edges(1:end-1) + diff(edges)/2;
            end
            pdf = ksdensity(values, centers, param{:});
            
        end
end

if flipbinvec
    centers = centers';
    edges = edges';
end
if flipvalues
    pdf = pdf';
end