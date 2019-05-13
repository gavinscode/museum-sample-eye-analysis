function [ mapResults, numReg ] = calculateRegionsInMap( vertices, indsIn, distThresh )
    %Written by Gavin Taylor, 2017. MIT License
    
    [  regionIds, ~, regVals]  = findRegionsOnMesh( vertices, indsIn, distThresh); 
    regionsToGet = find(regVals);
    mapResults = cell(length(regionsToGet),3);
    toKeep = [];
    numInRegions = [];
    
    for i = 1:length(regionsToGet)
        tempInds = zeros(length(indsIn),1);
        tempInds(regionIds == regionsToGet(i)) = 1;

        if sum(tempInds) > 3 %cut of is arbitary, but needs to be at least 3
            [lineVerticesT, borderIndsT, fullIndsT] = plotSortedIsoLine(vertices, tempInds, 0.5,0,distThresh);
            [lineVerticesT, sIndsT] = sortPointsWithTSP(lineVerticesT);
            borderIndsT = borderIndsT(sIndsT);
            lineVerticesT = vertcat(lineVerticesT,lineVerticesT(1,:));
            
            mapResults{i,1} = lineVerticesT; %close line...
            mapResults{i,2} = vertcat(borderIndsT,borderIndsT(1));
            mapResults{i,3} = fullIndsT;

            toKeep = [toKeep i];
            numInRegions = [numInRegions sum(tempInds)];
        end
    end
    
    mapResults = mapResults(toKeep,:);
    [~, sI] = sort(numInRegions,'descend');
    mapResults = mapResults(sI,:);
    
    numReg = size(mapResults,1);
end

