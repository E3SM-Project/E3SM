% edges for 3 hexagons

blue = [28,32,102]/255;
gray = [179,179,179]/255;
red = [163,0,0]/255;

[centers, edges, normals] = getHexCentersEdgesNormals([0,0], 1, 2);

plot(centers(:,1),centers(:,2),'.', 'MarkerSize', 25, 'MarkerEdgeColor', blue, 'MarkerFaceColor', blue);

plotEdges = zeros(size(edges,1)*3,2);
plotEdges(1:3:end,:) = squeeze(edges(:,1,:));
plotEdges(2:3:end,:) = squeeze(edges(:,2,:));
plotEdges(3:3:end,:) = NaN;

edgeCenters = 0.5*(squeeze(edges(:,1,:))+squeeze(edges(:,2,:)));

scale = 0.4*sign(normals(:,1)+1e-10*normals(:,2));
normals(:,1) = normals(:,1).*scale;
normals(:,2) = normals(:,2).*scale;

hold on;
plot(plotEdges(:,1),plotEdges(:,2), '-', 'Color', blue, 'LineWidth', 3);
%plot(edgeCenters(:,1), edgeCenters(:,2), '.', 'MarkerEdgeColor', gray, 'MarkerSize', 12);
blackMask = (edgeCenters(:,1).^2 + edgeCenters(:,2).^2) < 1;
for(index = find(~blackMask))
    start = [edgeCenters(index,1),edgeCenters(index,2)];
    stop = start+[normals(index,1),normals(index,2)];
    arrow(start,stop, 'EdgeColor', gray, 'FaceColor', gray, 'Width', 1.5, 'Length', 20);
end
for(index = find(blackMask))
    start = [edgeCenters(index,1),edgeCenters(index,2)];
    stop = start+[normals(index,1),normals(index,2)];
    arrow(start,stop, 'EdgeColor', 'k', 'FaceColor', 'k', 'Width', 1.5, 'Length', 20);
end
index = 17;
start = [edgeCenters(index,1),edgeCenters(index,2)];
stop = start+1.6*[normals(index,1),normals(index,2)];
plot(start(1),start(2),'.', 'MarkerSize', 25, 'MarkerEdgeColor', red, 'MarkerFaceColor', red);
arrow(start,stop, 'EdgeColor', red, 'FaceColor', red, 'Width', 1.5, 'Length', 15);

plot(.25,.25,'.', 'MarkerSize', 25, 'MarkerEdgeColor', red, 'MarkerFaceColor', red);
arrow([.25,.25],[.65,.3], 'EdgeColor', red, 'FaceColor', red, 'Width', 1.5, 'Length', 15);
hold off;
axis([-3, 3, -3, 3]); axis equal; axis off;

text('units','inch', 'position',[2.6 3.5], ...
    'fontsize',14, 'interpreter','latex', 'string',...
    '$$\hat{n}_i$$'); 

text('units','inch', 'position',[2.7 3.1], ...
    'fontsize',14, 'interpreter','latex', 'string',...
    '$$u_i$$'); 

text('units','inch', 'position',[3 2.85], ...
    'fontsize',14, 'interpreter','latex', 'string',...
    '$${\bf x}_j$$'); 

text('units','inch', 'position',[2.1 3.6], ...
    'fontsize',14, 'interpreter','latex', 'string',...
    '$${\bf x}_i$$'); 

text('units','inch', 'position',[3.3 3.15], ...
    'fontsize',14, 'interpreter','latex', 'string',...
    '$${\bf x}$$'); 

text('units','inch', 'position',[3.3 3.55], ...
    'fontsize',14, 'interpreter','latex', 'string',...
    '$${\bf u}_j({\bf x})$$'); 