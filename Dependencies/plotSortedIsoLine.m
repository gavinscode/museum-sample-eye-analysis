function [ lineVertices, lineIndices, highInds ] = plotSortedIsoLine(inputVertices, inputIntensity, contourThresh, sorting, localDist)
    %Written by Gavin Taylor, 2017. MIT License
    
%from main
highInds = find(inputIntensity >= contourThresh); 
lowInds = find(inputIntensity < contourThresh);
borderIndIt = 1;

if ~isempty(highInds) & ~isempty(lowInds)
    %work with the smallest list 
    %this loop is very slow for large arrays...
    if length(highInds) < length(lowInds)
        borderInds = zeros(length(highInds),1);
        for i = 1:length(highInds)
            dists = sqrt((inputVertices(:,1)-inputVertices(highInds(i),1)).^2+(inputVertices(:,2)-inputVertices(highInds(i),2)).^2+(inputVertices(:,3)-inputVertices(highInds(i),3)).^2);

            if ~isempty(localDist) %speeds up if provided
                tI = find(dists < localDist);
                [~, distInds] = sort(dists(tI));
                distInds = tI(distInds);
            else
                [~, distInds] = sort(dists);
            end
            % if nearest neihbros not all in high inds, then this is border
            if length(intersect(highInds, distInds(2:6))) ~= 5
                borderInds(borderIndIt) = highInds(i);
                borderIndIt = borderIndIt + 1;
            end
        end
        borderInds = borderInds(1:borderIndIt-1);
    else
        borderInds = zeros(length(lowInds),1);
        for i = 1:length(lowInds)
            dists = sqrt((inputVertices(:,1)-inputVertices(lowInds(i),1)).^2+(inputVertices(:,2)-inputVertices(lowInds(i),2)).^2+(inputVertices(:,3)-inputVertices(lowInds(i),3)).^2);

            if ~isempty(localDist)
                tI = find(dists < localDist);
                [~, distInds] = sort(dists(tI));
                distInds = tI(distInds);
            else
                [~, distInds] = sort(dists);
            end
            % if nearest neihbros not all in high inds, then this is border
            if length(intersect(lowInds, distInds(2:6))) ~= 5
                borderInds(borderIndIt) = lowInds(i);
                borderIndIt = borderIndIt + 1;
            end
        end
        borderInds = borderInds(1:borderIndIt-1);
    end

    % figure; plot3(inputVertices(highInds,1), inputVertices(highInds,2), inputVertices(highInds,3), 'or');

    %sorting usually now done externally with tsp
    if sorting
        %align to pca for finding angles
        tempVertices = [inputVertices(borderInds,1) - mean(inputVertices(borderInds,1)), inputVertices(borderInds,2) - mean(inputVertices(borderInds,2)), inputVertices(borderInds,3) - mean(inputVertices(borderInds,3))];
        pcaVecs = pca(tempVertices);
        rerToAxis = vectorRotationFromAtoB([0 0 1]',pcaVecs(:,3)');
        tempVertices = tempVertices*rerToAxis;
        angles = atan2(tempVertices(:,1), tempVertices(:,2));
        [~,  borderIndsSort] = sort(angles);
        borderIndsSort = borderInds(borderIndsSort);
        borderIndsSort = [borderIndsSort' borderIndsSort(1)];

        % figure; plot3(inputVertices(borderIndsSort,1), inputVertices(borderIndsSort,2), inputVertices(borderIndsSort,3));

        % %sort by smallest angular change about a certain point, also has skip over effects - is a poor solution to travelling salesman problem
        % centralPoint = mean(inputVertices(borderInds,:),1);
        % borderIndsSort = zeros(length(borderInds),1);
        % borderIndsSort(1) = borderInds(1);
        % for i = 2:length(borderInds);
        %     unusedInds = find(~isnan(borderInds));
        %     angles = zeros(length(unusedInds),1);
        %     vecA = inputVertices(borderIndsSort(i-1),:) - centralPoint;
        %     for j = 1:length(unusedInds)
        %         vecB = inputVertices(borderInds(unusedInds(j)),:) - centralPoint;
        %         angles(j) = atan2(norm(cross(vecA,vecB)),dot(vecA,vecB));
        %     end
        %     
        %     [~, nextInd] = min(angles);
        %     meanAngles(i) = sqrt((inputVertices(borderIndsSort(i-1),1)-inputVertices(borderInds(unusedInds(nextInd)),1)).^2+(inputVertices(borderIndsSort(i-1),2)-inputVertices(borderInds(unusedInds(nextInd)),2)).^2+...
        %         (inputVertices(borderIndsSort(i-1),3)-inputVertices(borderInds(unusedInds(nextInd)),3)).^2);
        %     borderIndsSort(i) = borderInds(unusedInds(nextInd));
        %     borderInds(unusedInds(nextInd))= NaN;
        % end
        %distance based sorting can cause  skip over effects
        % [~, maxDistInd] = max(sqrt((inputVertices(:,1)-mean(inputVertices(highInds(:),1))).^2+(inputVertices(:,2)-mean(inputVertices(highInds(:),2))).^2+(inputVertices(:,3)-mean(inputVertices(highInds(:),3))).^2));
        % borderIndsSort = zeros(length(borderInds),1);
        % borderIndsSort(1) = borderInds(1);
        % borderInds(1) = maxDistInd;
        % % figure; hold on
        % for i = 2:length(borderInds);
        %     dists = sqrt((inputVertices(borderIndsSort(i-1),1)-inputVertices(borderInds(:),1)).^2+(inputVertices(borderIndsSort(i-1),2)-inputVertices(borderInds(:),2)).^2+...
        %         (inputVertices(borderIndsSort(i-1),3)-inputVertices(borderInds(:),3)).^2);
        %     [~, nextInd] = min(dists);
        %     
        % %     plot3([headSphere(borderIndsSort(i-1),1) headSphere(borderInds(nextInd),1)], [headSphere(borderIndsSort(i-1),2) headSphere(borderInds(nextInd),2)],...
        % %         [headSphere(borderIndsSort(i-1),3) headSphere(borderInds(nextInd),3)], 'k-');
        % %     pause(0.1);
        %     borderIndsSort(i) = borderInds(nextInd);
        %     borderInds(nextInd) = maxDistInd;
        % end

        %solution to werid points which sometimes works
        % inds2Remove = find(meanAngles > mean(meanAngles)*1.1);
        % borderIndsSort(inds2Remove) = [];
        % borderIndsSort(end-2:end) = [];

        lineVertices = inputVertices(borderIndsSort,:);
        lineIndices = borderIndsSort;
    else
        lineVertices = inputVertices(borderInds,:);
        lineIndices = borderInds;
    end
else
    lineVertices = [];
    lineIndices = [];
    highInds = [];
end

end

