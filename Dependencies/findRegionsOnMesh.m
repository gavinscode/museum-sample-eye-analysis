function [ regionIds, regionlength, regionVals ] = findRegionsOnMesh( inputVertices, vertVals, localDist)
    %Written by Gavin Taylor, 2017. MIT License
    
    %similar to plot sorted isoline, but will assign points to regions
    %given (possibly similar) vertVals that are not connected

    if size(inputVertices,1) ~= size(vertVals,1)
       error('need vertVals and inputVertices to have same number of lines'); 
    end
    
    regionVals = [];
    regionlength = [];
    regionIndex = 0;
    regionIds = zeros(length(vertVals),1);
    indsAssigned = zeros(length(vertVals),1);
    indsToTest = zeros(length(vertVals),1);
    
    %intalise for first
    
%     plot3(inputVertices(:,1),inputVertices(:,2),inputVertices(:,3),'g*');
%     its = 1;
    %keep doing until all inds test
    while sum(indsAssigned) ~= length(vertVals)
       
       %if some inds remaining to test
       if sum(indsToTest)
           currInd = find(indsToTest);
           currInd = currInd(1);
           indsToTest(currInd) = 0;
           
%            if currInd == 576
%                 plot3(inputVertices(currInd,1),inputVertices(currInd,2),inputVertices(currInd,3),'g*');
%                 b = 1;
%            end
           
           if vertVals(currInd) == regionVals(regionIndex)
               dists = sqrt((inputVertices(:,1)-inputVertices(currInd,1)).^2+(inputVertices(:,2)-inputVertices(currInd,2)).^2+(inputVertices(:,3)-inputVertices(currInd,3)).^2);
               
               if ~isempty(localDist) %speeds up if provided
                    tI = find(dists < localDist);
                    [~, distInds] = sort(dists(tI));
                    distInds = tI(distInds);
               else
                    [~, distInds] = sort(dists);
               end
            
               regionlength(end) = regionlength(end) + 1;
               regionIds(currInd) = regionIndex;
               indsAssigned(currInd) = 1;
               for i = 2:6 %length(distInds)
                   if distInds(i) > length(indsAssigned) 
                       b = 1;
                   end
                  if indsAssigned(distInds(i)) == 0
                     indsToTest(distInds(i)) = 1;
                  end
               end
           end

       else
          %start new... 
          regionIndex = regionIndex + 1;
          
          currInd = find(indsAssigned == 0);
          currInd = currInd(1);
          
%           if abs(inputVertices(currInd,1) - 0.2641) < 0.05 & abs(inputVertices(currInd,2) - -0.3013) < 0.05 & abs(inputVertices(currInd,3) - -0.9162) < 0.05
%                b = 1;
%            end
          
          dists = sqrt((inputVertices(:,1)-inputVertices(currInd,1)).^2+(inputVertices(:,2)-inputVertices(currInd,2)).^2+(inputVertices(:,3)-inputVertices(currInd,3)).^2);
          if ~isempty(localDist) 
             tI = find(dists < localDist);
             [~, distInds] = sort(dists(tI));
             distInds = tI(distInds);
          else
             [~, distInds] = sort(dists);
          end
               
          regionlength(end + 1) = 1;
          regionIds(currInd) = regionIndex;
          regionVals(end + 1) = vertVals(currInd);
          indsAssigned(currInd) = 1;
          for i = 2:6;
              if indsAssigned(distInds(i)) == 0
                 indsToTest(distInds(i)) = 1;
              end
          end
       end
       
%        if indsToTest(576) == 1
%              plot3(inputVertices(currInd,1),inputVertices(currInd,2),inputVertices(currInd,3),'r*');
%              b = 1;
%        end
%        its = its + 1;
    end
end

