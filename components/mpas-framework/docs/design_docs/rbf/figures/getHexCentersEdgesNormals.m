function [cellCenters, cellEdges, cellNormals] = getHexCentersEdgesNormals(cellCenters, radius, maxDistance)

    function [edges] = getHexEdges(center)
        edges = zeros(6,2,2); % 6 edges, 2 ends to each edge, 2 dimensions
        theta = 2*pi*(0:6)/6;
        for(index = 1:6)
            edges(index,1,1) = center(1) + radius*cos(theta(index));
            edges(index,1,2) = center(2) + radius*sin(theta(index));
            edges(index,2,1) = center(1) + radius*cos(theta(index+1));
            edges(index,2,2) = center(2) + radius*sin(theta(index+1));
        end
    end

    function [centers] = getHexNieghbors(center)
        centers = zeros(6,2);
        theta = 2*pi*(0.5+(0:5))/6;
        
        for(index = 1:6)
            centers(index,1) = center(1) + sqrt(3)*radius*cos(theta(index));
            centers(index,2) = center(2) + sqrt(3)*radius*sin(theta(index));
        end
    end

newCenters = cellCenters;
cellEdges = [];

while(~isempty(newCenters))
    centersToAdd = [];
    for(centerIndex = 1:size(newCenters,1))
        newEdges = getHexEdges(newCenters(centerIndex,:));
        edgesToAdd = [];
        for(edgeIndex = 1:size(newEdges,1))
            if(isempty(cellEdges))
                edgesToAdd = newEdges(edgeIndex,:,:);
                continue;
            end
            distancesSquared = (cellEdges(:,1,1)-newEdges(edgeIndex,2,1)).^2 ...
                + (cellEdges(:,1,2)-newEdges(edgeIndex,2,2)).^2 ...
                + (cellEdges(:,2,1)-newEdges(edgeIndex,1,1)).^2 ...
                + (cellEdges(:,2,2)-newEdges(edgeIndex,1,2)).^2;
            if(min(distancesSquared) > 1e-5)
                edgesToAdd = [edgesToAdd; newEdges(edgeIndex,:,:)];
            end
        end
        cellEdges = [cellEdges; edgesToAdd];
        neighborCenters = getHexNieghbors(newCenters(centerIndex,:));
        for(neighborIndex = 1:length(neighborCenters))
            if((abs(neighborCenters(neighborIndex,1)) > maxDistance ) ...
                    || (abs(neighborCenters(neighborIndex,2)) > maxDistance ))
                continue;
            end
            distancesSquared = (cellCenters(:,1)-neighborCenters(neighborIndex,1)).^2 ...
                + (cellCenters(:,2)-neighborCenters(neighborIndex,2)).^2;
            if(min(distancesSquared) > 1e-5)
                centersToAdd = [centersToAdd; neighborCenters(neighborIndex,:)];
            end
        end
    end
    newCenters = centersToAdd;
    cellCenters = [cellCenters; centersToAdd];
end

cellNormals = zeros(size(cellEdges,1),2);

cellNormals(:,1) = cellEdges(:,2,2)-cellEdges(:,1,2);
cellNormals(:,2) = -(cellEdges(:,2,1)-cellEdges(:,1,1));
scale = 1./sqrt(cellNormals(:,1).^2 + cellNormals(:,2).^2);
cellNormals(:,1) = cellNormals(:,1).*scale;
cellNormals(:,2) = cellNormals(:,2).*scale;

% figure(1);
% hold off;
% for(index = 1:size(cellEdges,1))
%     plot(squeeze(cellEdges(index,:,1)),squeeze(cellEdges(index,:,2)), '-', 'LineWidth', 2, 'Color', [15,111,198]/255);
%     hold on;
% end
% plot(cellCenters(:,1),cellCenters(:,2), '.', 'MarkerSize', 15, 'MarkerEdgeColor', [0.4,0.4,0.4]);
% hold off;
% axis([-3,3,-3,3]); axis equal;
% axis off;

end