function [ norm ] = normalToSurface( pointCoords, surface, referenceCenter, flipAway, rng)
    %Written by Gavin Taylor, 2017. MIT License

    %note this now forces an even poly and centers x...
    
    %flip away should make normal point into a convex region, if zero it
    %should point away from the convex region...

    %_withCurvature fiel implements flat feature checking and also forces
    %changes on border
    
    %mostly pulled from AngleOfRefractedRay

    debugPlot = 0;
    
%     rng = 50;    %this may need tuning...
    if ~isempty(rng)
        surfInds = find(sqrt((surface(:,1)-pointCoords(1)).^2 + (surface(:,2)-pointCoords(2)).^2 + (surface(:,3)-pointCoords(3)).^2) < rng);
    else
       %use all inds
       surfInds = 1:size(surface,1);
    end
    
    if debugPlot
         figure; subplot(1,3,1); hold on; axis equal
         plot3(surface(surfInds,1), surface(surfInds,2), surface(surfInds,3), '.');
        
         line([0 0], [0 10], [0 0],'color', 'r');
         line([0 10], [0 0], [ 0 0],'color', 'g');
         line([0 0], [0 0], [0 10],'color', 'b');
    end
    
    referenceCenterOrig = referenceCenter;
    
    %rotate surface to axis of max variation
    centM = mean(surface(surfInds,:));
    referenceCenter = referenceCenter-centM;
    
    surfPts = surface(surfInds,:)-centM;
%     surfPts = [surface(surfInds,1)-pointCoords(1),surface(surfInds,2)-pointCoords(2),surface(surfInds,3)-pointCoords(3)];
    pcaVecs = pca(surfPts);
    
    rer = vectorRotationFromAtoB([0 0 1]',pcaVecs(:,3)');
    referenceCenter = referenceCenter*rer;
    %check rotation is in correct direction
    if referenceCenter(3) < 0
        rer = vectorRotationFromAtoB([0 0 -1]',pcaVecs(:,3)');
        rer2 = vectorRotationFromAtoB(pcaVecs(:,3)',[0 0 -1]');
        referenceCenter = referenceCenterOrig-centM;
        referenceCenter = referenceCenter*rer;
    else
        rer2 = vectorRotationFromAtoB(pcaVecs(:,3)',[0 0 1]');  
    end
    
    newPts = surfPts*rer;
    
    %shift base to origin - required, as intersect should already be close to zero?
    [~, centInd] = min(sqrt(newPts(:,1).^2+newPts(:,1).^2)); 
    newPts(:,3) = newPts(:,3) - newPts(centInd,3);
    
    if debugPlot
         subplot(1,3,2); hold on; axis equal
         plot3(newPts(:,1), newPts(:,2), newPts(:,3), '.');
         
         xVec = [10 0 0]*rer;
         yVec = [0 10 0]*rer;
         zVec = [0 0 10]*rer;
         
         line([0 xVec(1)], [0 xVec(2)], [0 xVec(3)],'color', 'r');
         line([0 yVec(1)], [0 yVec(2)], [0 yVec(3)], 'color', 'g');
         line([0 zVec(1)], [0 zVec(2)], [0 zVec(3)],'color', 'b');
    end
    
    %fit the plane
    %force to touch origin
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    %forced to touch zero
%     opts.Lower = [0 -Inf -Inf -Inf -Inf -Inf]; %
%     opts.Upper = [0 Inf Inf Inf Inf Inf]; 
    %for even
%     opts.Lower = [-Inf 0 0 -Inf 0 -Inf]; %
%     opts.Upper = [Inf 0 0 Inf 0 Inf];

    %%% currently used!
    %   for force concave around lens - x*y can be weird if limited area
    %   but then need to make sure orientation is correct
    opts.Lower = [-Inf -Inf -Inf 0 0 0]; %
    opts.Upper = [Inf Inf Inf Inf Inf Inf];

    fittedPoly = fit(newPts(:,1:2), newPts(:,3),'poly22', opts);      %very slow - could pass it a downsampled array?
    
    if debugPlot
        plot(fittedPoly);
    end
    
    %find the closest point to zero
    function F = solveTheFn(in,poly)
        F = sqrt(poly(in(1),in(2))^2+in(1)^2+in(2)^2);
    end
    
    %fit at zero should be best if even polynomial used
    closeP = [0 0];

    [dX, dY] = differentiate(fittedPoly, closeP(1), closeP(2));
    if isnan(dX)
        dX = 0;
    end
    if isnan(dY)
        dY = 0;
    end
    N = [dX, dY, -1]; N = N/sqrt(N(1)^2+N(2)^2+N(3)^2);
    
    %flip normal to point towards lens
    referenceCenter = [referenceCenter(1) - pointCoords(1), referenceCenter(2) - pointCoords(2), referenceCenter(3) - pointCoords(3)];
    referenceCenter = referenceCenter*rer;
    referenceCenter = referenceCenter/sqrt(referenceCenter(1)^2+referenceCenter(2)^2+referenceCenter(3)^2);
    if ~flipAway
        if sqrt((referenceCenter(1)-N(1))^2+(referenceCenter(2)-N(2))^2+(referenceCenter(3)-N(3))^2) > sqrt(2)
            N = N*-1;
        end
    else
        if sqrt((referenceCenter(1)-N(1))^2+(referenceCenter(2)-N(2))^2+(referenceCenter(3)-N(3))^2) < sqrt(2)
            N = N*-1;
        end
    end
    
    if debugPlot
        line([closeP(1) N(1)*50],[closeP(2) N(2)*50] , [0 N(3)*50]);
    end
    
    norm = N*rer2;
    
     if debugPlot
         subplot(1,3,1);
         line([pointCoords(1) pointCoords(1)+norm(1)*50],[pointCoords(2) pointCoords(2)+norm(2)*50] , [pointCoords(3) pointCoords(3)+norm(3)*50]);
    end
end
