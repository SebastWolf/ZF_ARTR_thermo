function [dist,path,pred] = shortest_path_free_energy(startloc, endloc, free_energy_all)


regionsize = size(free_energy_all);

R = regionsize(1);
C = regionsize(2);

startloc = startloc(:)';
endloc = endloc(:)';

% create a graph that connects all nodes to their immediate neighbors -
% left, right, up or down.
[y,x] = meshgrid(1:C,1:R);
xy0 = [x(:),y(:)];
xy1 = [xy0 + [0 1];xy0 + [1 0];xy0 - [0 1];xy0 - [1 0]];
xy0 = repmat(xy0,4,1);

% discard edges that wander out of the region boundary
dropind = (xy1(:,1) < 1) | (xy1(:,2) < 1) | (xy1(:,1) > R) | (xy1(:,2) > C);
xy1(dropind,:) = [];
xy0(dropind,:) = [];

% list of all valid edges in the graph
edges = [sub2ind(regionsize,xy0(:,1),xy0(:,2)),...
    sub2ind(regionsize,xy1(:,1),xy1(:,2))];

% convert startloc and endloc into nodes
% since the edge list contains all possible edges in the rectangular grid,
% as long as the start and end nodes both live inside the grid, then there
% will always be at least one path between them.
startnode = sub2ind(regionsize,startloc(1),startloc(2));
endnode = sub2ind(regionsize,endloc(1),endloc(2));

% define the weigths matrix
clear W
for i = 1:size(edges,1);
    W(i,1) = free_energy_all(edges(i,2))-free_energy_all(edges(i,1));
end
%W(W>0) = 1e40;
W = abs(W);

% Create a directed graph with 6 nodes and 11 edges.
%W = [.41 .99 .51 .32 .15 .45 .38 .32 .36 .29 .21];
DG = sparse(edges(:,1),edges(:,2),W);
%h = view(biograph(DG,[],'ShowWeights','on'))


% Find the shortest path
[dist,path,pred] = graphshortestpath(DG,startnode,endnode,'Method','Bellman-Ford');

% Mark the nodes and edges of the shortest path by coloring them red and
% increasing the line width.
% set(h.Nodes(path),'Color',[1 0.4 0.4])
% edges = getedgesbynodeid(h,get(h.Nodes(path),'ID'));
% set(edges,'LineColor',[1 0 0])
% set(edges,'LineWidth',1.5)
end
