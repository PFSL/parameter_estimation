function plot_graph(nodes, elem, ngl)

%% simple routine to plot lines describing elements of truss

[nelem,~] = size(elem);     % number of elements

% set limits of plot for simple truss example
xmin = min(nodes(:,2));
xmax = max(nodes(:,2));
ymin = min(nodes(:,3));
ymax = max(nodes(:,3));
set(gca, 'XLim', [xmin - 1, xmax + 1], 'YLim', [ymin - 1, ymax + 1]);

figure(1);  % define a new figure
% loop over elements to plot using line command
for i=1:nelem
    
    % nodes
    noi = elem(i,2);
    nof = elem(i,3);
    % initial coordinates
    xi = nodes(noi,2);
    yi = nodes(noi,3);
    % final coordinates
    xf = nodes(nof,2);
    yf = nodes(nof,3);
    % coordinate vectors
    x = [[xi],[xf]];
    y = [[yi],[yf]];
    
    s = line(x', y', 'Color', 'b', 'LineWidth', 1);
    
end

end