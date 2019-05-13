%Written by Gavin Taylor, 2017. MIT License

% clear
clc ;close all

%%% flow can be controlled can be switched using 'do' flags for loading,
%%% calculating and plotting (and saving data or images) in configeration
%%% section


%--- configeration
doLoading = 1;
    %%% this is configured to use data from the Temperate gnatt folder as an example

    fileFolder = 'Enter directory to load data from';
%best to load stl files into meshlab to clean up. Save as little endian stl files from Amira. Do not unify duplicated vertices in meshlab.
    eyeSTLFile = 'Te_LeftEye.stl';
    headSTLFile = []; %optional, just for plotting
%both lines and paths come from an amira linset file. Edit in VIM to save vertices as lines and index as path (seperate files)
    eyeLineFile = {'Te_LeftEyeLines.am'}; %if linesets become split in Amira, enter the name for each as a seperate cell
    eyePathfile = {'Te_LeftEyePaths.am'};
%transformations
    %get 4x4 transfrom using getTransform command in amira
    eyeTransform = [0.571225 0.637066 -0.517543 0 0.797134 -0.280285 0.534809 0 0.195649 -0.718047 -0.667932 0 -533.632 428.023 329.986 1]; 
    headTransform = [0.571225 0.637066 -0.517543 0 0.797134 -0.280285 0.534809 0 0.195649 -0.718047 -0.667932 0 -533.632 428.023 329.986 1];
    rightEyeTransform = [0.994715 -0.102642 0.00295139 0 0.102544 0.991444 -0.0807688 0 0.00536428 0.0806446 0.996729 0 -20.7587 25.7959 12.8311 1];
    
    %get volume sizes from crop editor in amira (after resample transformed)
    rightEyeVolSize = [508 491 513];
    %get min coord from crop editor
    rightEyeMinCoord = [11.2981 -2.65581 -203.681];
    
    %voxel size of the orignal volume
    voxelSize = 1; %of original volumes, should be isotropic and same for all
    autoAlignPitchBetweenEyes = 1;
    %%% could also add auto-alignment code to get head square in world axis, but not sure this is meaningfull given different eye shapes.

doCalculations = 1;
    smoothingFactor = []; %weight of smoothing for local vectors and lens diameters, empty will not smooth
    
    %%%re-think this limit
    stdLimitForAcceptanceAng = 4; %constrains unrealistically large acceptance angles
    %%%choose acceptance angle such that V-cutoff as close as possible to Vs for given intensity*contrast
    
    sphereDensity = 5; %num points on test sphere (3->642, 4->2562, 5->10242, 6->40962,...
    imageHeight = 500; %for drawing in unwrapped images (in pixels), width will be double
    useLargestRegion = 0; %only plot largest indivdual regions of fov
    
doPlotting = 1;
    plotSphereR = 10^2.8; %set such that it is a bit bigger than the head..
    bgColor = 1; %0 black or 1 white
    numColBins = []; %if empty, try to use preloaded values from MergeAndPlot
    plotBinoLine = 1;
    plotRightEyeLine = 0; %will also plot left line in some cases
    displayTitles = 0;
    displayAxisLabelsOnUnWrappedImages = 1;
    mSize =[]; %marker size of bino facets on eye - set to zero to turn off
    
    %to display contour enter: 1 - Lens Diameter, 2 - Interommatidial angle, 3 - Lens parameter
    %displaying more than one contour per plot will probably be hard to read
    %%% hard to read with more than ~3 contours
    dispContourOnDiameterImage = [2];
    dispContourOnInterOImage = [1];
    dispContourOnLensPImage = [2];
    
    %must have run histogram analysis on multiple eyes to generate this
    useGroupColBins = 1;
    
    simUserImage = 1; %otherwise will use checks
    %%%% note: there is no option to display color images, but this is easy to implement
    imageToGet = 'shutterstock_651336862.jpg';
    imType = 'JPEG';
    imgRotate = 0; %positive is left, negative right
    contrastAdj = 0.2; %scales image from [x 1-x]

    numChequers = 20; %#checks

    sensitivityScale = 0; %just use 1 of these
    sensitvityScaleByColSat = 0;
        invertSensColScale = 0;
    
    autoAlignAndIllum = 1;
        AlignAngle = [90 -30]; %front on
        
    simBothEyes = 1; %otherwise just plot left

    saveFileName = 'EyeData';
    
saveData = 1;
    saveDirectory = 'Enter directory to save data';
%%% When data is saved use MergeAndPlot.m to plot combined histograms 
% and also to generate color bins across variable ranges for multiple eyes
    
saveImages = 1;
    saveImDirectory = 'Enter directory to save images';

%--- loading
if doLoading
    cd(fileFolder)
    
    %load surfaces
    [F,V,eyeNorms] = stlread(eyeSTLFile); 
    if ~isempty(headSTLFile)
        [FHead,VHead,HeadNorms] = stlread(headSTLFile); 
    end
    
    %uniquify verts
    [V, cutInds, fullInds] = unique(V, 'rows');
    F = fullInds(F);
    
    % load lines and paths
    if length(eyeLineFile) == 1
        %get single path
        line1 = {load(eyeLineFile{1})};
        path1 = {load(eyePathfile{1})};
    else
        %get multiple paths :(
        line1 = cell(length(eyeLineFile),1);
        path1 = cell(length(eyeLineFile),1);
        for i = 1:length(eyePathfile)
            line1{i} = load(eyeLineFile{i});
            path1{i} = load(eyePathfile{i});
        end
    end
    
    %shift up all inds (must only do once)
    totalPaths = 0;
    for i = 1:length(eyePathfile)
        pathT = path1{i};
        if sum(pathT == -1) > 0 
            pathT = pathT + 1;
            pathT(isnan(pathT)) = [];
        end
        totalPaths = totalPaths + sum(pathT == 0);
        path1{i} = pathT;
    end
    
    %parse paths
    eachPaths = cell(totalPaths,1);
    
    iterator = 1;
    for i = 1:length(eyePathfile)
        pathT = path1{i};
        lineT = line1{i};
        changeInds = [0 find(pathT == 0)'];
        for j = 1:length(changeInds)-1
            inds = unique(pathT(changeInds(j)+1:changeInds(j+1)-1));
            eachPaths{iterator} = lineT(inds,:);
            iterator = iterator + 1;
        end
    end
    
    %
    headTransform = reshape(headTransform, [4,4]);
    eyeTransform = reshape(eyeTransform, [4,4]);
    rightEyeTransform = reshape(rightEyeTransform, [4,4]);
    
    %new method of adjusting vertices
    if ~isempty(headSTLFile)
        VHeadShifted = [VHead ones(size(VHead(:,3),1),1)]*headTransform';
        VHeadShifted = VHeadShifted(:,1:3);
    end
        
    Vshifted = [V ones(size(V(:,3),1),1)]*eyeTransform';
    Vshifted = Vshifted(:,1:3);
    
    VshiftedMirror = Vshifted-rightEyeMinCoord;
    VshiftedMirror(:,1) = -VshiftedMirror(:,1)+(rightEyeVolSize(1)/voxelSize-1)*voxelSize;
    VshiftedMirror = VshiftedMirror+rightEyeMinCoord;
    VshiftedMirror = [VshiftedMirror ones(size(VshiftedMirror(:,3),1),1)]*rightEyeTransform';
    VshiftedMirror = VshiftedMirror(:,1:3);

    meanPosLeft = mean(Vshifted); pcaLeft = pca(Vshifted-meanPosLeft);
    meanPosRight = mean(VshiftedMirror); pcaRight = pca(VshiftedMirror-meanPosRight);
    
    if autoAlignPitchBetweenEyes
        %adjust pitch of right eye to match left eye
        leftPitch = atan2(pcaLeft(2,2),pcaLeft(3,2));
        rightPitch = atan2(pcaRight(2,2),pcaRight(3,2));
        pitchRotAdj = vrrotvec2mat([1 0 0 leftPitch-rightPitch]);
        VshiftedMirror = (VshiftedMirror-meanPosRight)*pitchRotAdj+meanPosRight;
        
%         figure; hold on; axis equal
%         plot3(Vshifted(:,1), Vshifted(:,2), Vshifted(:,3));
%         plot3(VshiftedMirror(:,1), VshiftedMirror(:,2), VshiftedMirror(:,3));
%         line([0 pcaLeft(1,1)]*100, [0 pcaLeft(2,1)]*100, [0 pcaLeft(3,1)]*100, 'color','r','linewidth',2);
%         line([0 pcaLeft(1,2)]*100, [0 pcaLeft(2,2)]*100, [0 pcaLeft(3,2)]*100, 'color','r','linewidth',4);
%         line([0 pcaLeft(1,3)]*100, [0 pcaLeft(2,3)]*100, [0 pcaLeft(3,3)]*100, 'color','r','linewidth',1);
%         line([0 pcaRight(1,1)]*100, [0 pcaRight(2,1)]*100, [0 pcaRight(3,1)]*100, 'color','b','linewidth',2);
%         line([0 pcaRight(1,2)]*100, [0 pcaRight(2,2)]*100, [0 pcaRight(3,2)]*100, 'color','b','linewidth',4);
%         line([0 pcaRight(1,3)]*100, [0 pcaRight(2,3)]*100, [0 pcaRight(3,3)]*100, 'color','b','linewidth',1);
    end
    
    headCenter = mean([Vshifted', VshiftedMirror']');
    
    %shift to head center
    Vshifted = [Vshifted(:,1)-headCenter(1), Vshifted(:,2)-headCenter(2), Vshifted(:,3)-headCenter(3)];
    VshiftedMirror = [VshiftedMirror(:,1)-headCenter(1), VshiftedMirror(:,2)-headCenter(2), VshiftedMirror(:,3)-headCenter(3)];
    if ~isempty(headSTLFile)
        VHeadShifted = [VHeadShifted(:,1)-headCenter(1), VHeadShifted(:,2)-headCenter(2), VHeadShifted(:,3)-headCenter(3)];
    end
%     if autoAlignBothEyes
%         pcaAx = pca([Vshifted' VshiftedMirror']');
%     end
    
    %rotate normals to get nice visualization
    if ~isempty(headSTLFile)
        HeadNormsRot = [HeadNorms zeros(size(HeadNorms(:,3),1),1)]*headTransform';
        HeadNormsRot = HeadNormsRot(:,1:3);
    end
    
    eyeNormsShiftedRot = [eyeNorms zeros(size(eyeNorms(:,3),1),1)]*eyeTransform';
    eyeNormsShiftedRot = eyeNormsShiftedRot(:,1:3);

    figure; hold on; axis vis3d;
    patch(struct('faces', F,'vertices',Vshifted),'FaceColor', 'b', 'EdgeColor', 'none'); %[0 0 0]
    patch(struct('faces', F,'vertices',VshiftedMirror),'FaceColor', 'r', 'EdgeColor', 'none'); %[0 0 0]
    if ~isempty(headSTLFile)
        patch(struct('faces', FHead,'vertices',VHeadShifted),'FaceColor', 'g', 'EdgeColor', 'none');
    end
    xlabel('x'); ylabel('y'); zlabel('z');
    title('Check that eyes correctly located on head adjust transform if not');
end

%--- calculations
if doCalculations
    %assign patches and calculate vectors from eyes

    %values saved
    eachPatch = cell(totalPaths,2);
    patchTips = zeros(totalPaths,3,2);
    patchNorms = zeros(totalPaths,3,2);

    %working values when assigning patches
    vertLabels = zeros(size(V,1),1);
    ownRadius = zeros(totalPaths,1);

    for i = 1:totalPaths;%size(F,1)
        %center and take pca of path
        tempPath = eachPaths{i};
        meanPt = mean(tempPath);
        tempPath = [tempPath(:,1)-meanPt(1), tempPath(:,2)-meanPt(2), tempPath(:,3)-meanPt(3)];
        pcaMat = pca(tempPath); 
        rotMat = vectorRotationFromAtoB([0 0 1],pcaMat(:,3));

        %note that paths are fitted to orignal untransformed version of eye.
        %apply to all the original vertices in eye
        tempFulVerts = [V(:,1)-meanPt(1), V(:,2)-meanPt(2), V(:,3)-meanPt(3)];
        tempFulVerts = tempFulVerts*rotMat;
        tempPath = tempPath*rotMat;

        %find all the pts within x,y radius
        ownRadius(i) = mean(sqrt(tempPath(:,1).^2 + tempPath(:,2).^2));
        fullDists = sqrt(tempFulVerts(:,1).^2 + tempFulVerts(:,2).^2);

        %shouldnt be too far away in z either, or already labelled
        closeInds = find(fullDists < ownRadius(i) & abs(tempFulVerts(:,3)) < 50 & vertLabels == 0);

        %flag vertices found
        vertLabels(closeInds) = i;

        %copy patches from each set of aligned vertices
        eachPatch{i,1} = Vshifted(closeInds,:);
        eachPatch{i,2} = VshiftedMirror(closeInds,:); %indices on original eye correspond to mirrored eye

        %calculate facet surf center and pca, then shift and rotate to center point
        facetVerts = Vshifted(closeInds,:);
        meanPt = mean(facetVerts);
        tempFacetVerts = [facetVerts(:,1)-meanPt(1), facetVerts(:,2)-meanPt(2), facetVerts(:,3)-meanPt(3)];
        pcaMat = pca(tempFacetVerts); 
        rotMat = vectorRotationFromAtoB([0 0 1],pcaMat(:,3));
        tempFacetVerts = tempFacetVerts*rotMat;

        %record center indices and normal
        [~, centerInd] = min(sqrt(tempFacetVerts(:,1).^2 + tempFacetVerts(:,2).^2));
        patchTips(i,:,1) = facetVerts(centerInd,:);
        patchNorms(i,:,1) = normalToSurface( facetVerts(centerInd,:), facetVerts, meanPt, 1,[]);

        %repeat for mirrored vertices
        facetVerts = VshiftedMirror(closeInds,:);
        meanPt = mean(facetVerts);
        tempFacetVerts = [facetVerts(:,1)-meanPt(1), facetVerts(:,2)-meanPt(2), facetVerts(:,3)-meanPt(3)];
        pcaMat = pca(tempFacetVerts); 
        rotMat = vectorRotationFromAtoB([0 0 1],pcaMat(:,3));
        tempFacetVerts = tempFacetVerts*rotMat;
        [~, centerInd] = min(sqrt(tempFacetVerts(:,1).^2 + tempFacetVerts(:,2).^2));

        patchTips(i,:,2) = facetVerts(centerInd,:);
        patchNorms(i,:,2) = normalToSurface( facetVerts(centerInd,:), facetVerts, meanPt, 1,[]);
    end

    emptyInds = find(vertLabels == 0);

    %- link the patches into neighborhoods
    %linking is valid for left and right eyes

    patchConnections = cell(totalPaths,1);

    %find patches that have parts of their patch close to another one
    for i = 1:totalPaths
        dists = sort(sqrt((patchTips(i,1,1) - patchTips(:,1,1)).^2 +(patchTips(i,2,1) - patchTips(:,2,1)).^2 +(patchTips(i,3,1) - patchTips(:,3,1)).^2));
        dists(1) = []; %remove 1st as it is owner

        %uses distance to nearest neighbor and fraction of own radius as test for proximity - generally works
        proxInds = find(sqrt((Vshifted(:,1) - patchTips(i,1,1)).^2 + (Vshifted(:,2) - patchTips(i,2,1)).^2 + ...
            (Vshifted(:,3) - patchTips(i,3,1)).^2) < dists(1)+ownRadius(i)*0.25);

        %copy labels of proximal vertices
        nearbyLabels = vertLabels(proxInds);
        nearbyLabels(nearbyLabels == i) = [];
        nearbyLabels(nearbyLabels == 0) = [];

        %store as neighbors
        patchConnections{i} = unique(nearbyLabels);
    end

    %check each patch is connected to all of its connections
    for i = 1:totalPaths
        tempConns = patchConnections{i};
        for j = 1:length(tempConns)
            patchConnections{tempConns(j)} = unique([patchConnections{tempConns(j)}' i]');
        end
    end

    %check patches are not linked more than 6 times
    numConnections = zeros(totalPaths,1);
    for i = 1:totalPaths
        numConnections(i) = size(patchConnections{i},1);

        %remove furthest if more than 6
        if numConnections(i) > 6
            tempConns = patchConnections{i};

            %sort distance to center of all neighbors
            [~,SInd] = sort(sqrt((patchTips(i,1,1) - patchTips(patchConnections{i},1,1)).^2 + ...
                (patchTips(i,2,1) - patchTips(patchConnections{i},2,1)).^2 +...
                (patchTips(i,3,1) - patchTips(patchConnections{i},3,1)).^2));

            patchConnections{i} = tempConns(SInd(1:6));
            numConnections(i) = 6;

            %delete from former neighbors
            for j = 7:length(SInd)
                testCons = patchConnections{tempConns(SInd(j))};
                testCons(testCons == i) = [];
                patchConnections{tempConns(SInd(j))} = testCons;
                numConnections(tempConns(SInd(j))) = length(testCons);
            end
        end
    end
    
    %- calculate parameters from patches
    %calculations can be used for both eyes

    if ~isempty(smoothingFactor)
        %locally smooth vectors based on their neighbors
        patchNormsTemp = zeros(size(patchNorms));
        for i = 1:totalPaths
            tempConns = patchConnections{i};

            patchNormsTemp(i,:,1) = wmean(patchNorms([i tempConns'],:,1),[ones(1,3)' smoothingFactor*ones(3,length(tempConns))]');
            patchNormsTemp(i,:,2) = wmean(patchNorms([i tempConns'],:,2),[ones(1,3)' smoothingFactor*ones(3,length(tempConns))]');
        end
        patchNorms = patchNormsTemp;
    end

    % calculate angles & diameters
    interOAngle = zeros(totalPaths,1);
    lensDiamteres = zeros(totalPaths,1);

    for i = 1:totalPaths
        tempConns = patchConnections{i};

        %calcaulte and average interommatidial angles from neighbors normals
        angles = zeros(length(tempConns),1);
        for j = 1:length(tempConns)
            angles(j) = acos(dot(patchNorms(i,:,1),patchNorms(tempConns(j),:,1))/(norm(patchNorms(i,:,1))*norm(patchNorms(tempConns(j),:,1))));
        end
        interOAngle(i) = mean(angles)/pi*180;

        %diameter for CE facet analysis is average distance to neighbors 
        lensDiamteres(i) = mean(sqrt((patchTips(i,1,1)-patchTips(tempConns,1,1)).^2 + ...
            (patchTips(i,2,1)-patchTips(tempConns,2,1)).^2 + (patchTips(i,3,1)-patchTips(tempConns,3,1)).^2 ));

        %optical radius of eye its diameter divided by the interommatidial angle
    end

    if ~isempty(smoothingFactor)
        % smooth lens diameters - angles effectively averaged by smoothing normals
        lensDsTemp = zeros(size(lensDiamteres));
        for i = 1:totalPaths
            goodIds = patchConnections{i};
            lensDsTemp(i) = wmean(lensDiamteres([i goodIds']),[1 smoothingFactor*ones(1,length(goodIds))]');
        end
        lensDiamteres = lensDsTemp;
    end

    % calculate optical parameters for facets
    acceptanceAngs = zeros(length(interOAngle),1);
    sensitivites = zeros(length(interOAngle),1);

    % calculate acceptance angle
    for i = 1:length(interOAngle)
        acceptanceAngs(i) = interOAngle(i);
        %can also calcualte based on optimal and diffraction, but may be unrealistic
%         acceptanceAngs(i) = sqrt((1.06*interOAngle(i)/180*pi)^2+(0.500/lensDiamteres(i))^2)/pi*180;
    end

    % limit acceptance angles to be some multiple of their standard deviation from the mean
    %       on eyes with few facets, sometimes the interommatidial angle is very
    %       large on the periphery creating unrealistically large angles.
    cutAbove = mean(acceptanceAngs)+std(acceptanceAngs)*stdLimitForAcceptanceAng;
    acceptanceAngs(acceptanceAngs > cutAbove) = cutAbove;

    %based on facet size and assumed acceptance angle - rhabdom length would influence this
    sensitivites = lensDiamteres.^2.*(acceptanceAngs/180*pi).^2;
    
    % calculate senstivity 
%     for i = 1:length(acceptanceAngs)
%         %preferably use acceptance angle of rhabdom (calculated from interommatidial angle) without diffraction
%         if acceptanceAngs(i) < cutAbove
%             sensitivites(i) = lensDiamteres(i)^2*(1.06*interOAngle(i)/180*pi)^2;
%         else
%             %have to use entire thing if acceptance angle was limited
%             sensitivites(i) = lensDiamteres(i)^2*(acceptanceAngs(i)/180*pi)^2;
%         end
%     end
    
    eyeParam = (interOAngle/180*pi.*lensDiamteres);
        
    %- calculate intersections to world sphere and fields of view
    intersectionPoint = zeros(totalPaths,3,2);

    r = 10^9; %radius of test sphere, should be effectively infinity

    %step through for each facet
    for j = 1:totalPaths
        %calculate intersctiong of normal from facet to sphere of given radius
        origin = patchTips(j,:,1); 
        normer = patchNorms(j,:,1);
        p = (2*normer(1)* origin(1)+2*normer(2)*origin(2)+2*normer(3)*origin(3))/(normer(1)^2+normer(2)^2+normer(3)^2);
        q = (origin(1)^2+origin(2)^2+origin(3)^2-r^2)/(normer(1)^2+normer(2)^2+normer(3)^2);
        val = -(p/2)+sqrt((p/2)^2-q);
        lensNormsVal = normer.*val;
        intersectionPoint(j,:,1) = [lensNormsVal(1)+origin(1), lensNormsVal(2)+origin(2), lensNormsVal(3)+origin(3)];

        %do same for mirrored eye
        origin = patchTips(j,:,2); 
        normer = patchNorms(j,:,2);
        p = (2*normer(1)* origin(1)+2*normer(2)*origin(2)+2*normer(3)*origin(3))/(normer(1)^2+normer(2)^2+normer(3)^2);
        q = (origin(1)^2+origin(2)^2+origin(3)^2-r^2)/(normer(1)^2+normer(2)^2+normer(3)^2);
        val = -(p/2)+sqrt((p/2)^2-q);% +-?
        lensNormsVal = normer.*val;
        intersectionPoint(j,:,2) = [lensNormsVal(1)+origin(1), lensNormsVal(2)+origin(2), lensNormsVal(3)+origin(3)];
    end
    intersectionPoint = intersectionPoint/r;
    
    %create a spherical mesh to map points onto
    TR = SubdivideSphericalMesh(IcosahedronMesh,sphereDensity);
    headSphere = TR.X;

    %find general minimum distance between points
    minDst = 0;
    for i = round(1:size(headSphere,1)/50:size(headSphere,1))
        toTest = 1:size(headSphere,1); toTest(i) = [];
        minDst = minDst + min(sqrt((headSphere(i,1) - headSphere(toTest,1)).^2 + (headSphere(i,2) - headSphere(toTest,2)).^2 + (headSphere(i,3) - headSphere(toTest,3)).^2));
    end
    minDst = minDst/length(1:100:size(headSphere,1));

    %find closest point on mesh for each the intersections of each eye
    closestPoint = zeros(totalPaths,1);
    for i = 1:totalPaths
        dists = sqrt((intersectionPoint(i,1,1)-headSphere(:,1)).^2+(intersectionPoint(i,2,1)-headSphere(:,2)).^2+(intersectionPoint(i,3,1)-headSphere(:,3)).^2);
        [~, minInd] = min(dists);
        closestPoint(i) = minInd;
    end

    closestPointRight = zeros(totalPaths,1);
    for i = 1:totalPaths
        target = intersectionPoint(i,:,2);
        dists = sqrt((target(1)-headSphere(:,1)).^2+(target(2)-headSphere(:,2)).^2+(target(3)-headSphere(:,3)).^2);
        [~, minInd] = min(dists);
        closestPointRight(i) = minInd;
    end

    %calcualte borders of visual field
    %test what points on sphere are inside border of points
    faceSmooth = 0.9; %may need tuning, ideally 1
    closerPts = (intersectionPoint(:,:,1))*0.9;
    furtherPts = (intersectionPoint(:,:,1))*1.1;
    fullPts = [closerPts'  furtherPts']'; %create a 'box' around intersection inds
    bFace = boundary(fullPts(:,1), fullPts(:,2), fullPts(:,3),faceSmooth); %triangulate boundary of box
    inFOVLeft = double(intriangulation(fullPts, bFace, headSphere)); %calculate points on sphere in the triangulation

    %do other side
    closerPts = (intersectionPoint(:,:,2))*0.9;
    furtherPts = (intersectionPoint(:,:,2))*1.1;
    fullPts = [closerPts' furtherPts']'; 
    bFace = boundary(fullPts(:,1), fullPts(:,2), fullPts(:,3),faceSmooth);
    inFOVRight = double(intriangulation(fullPts, bFace, headSphere));

    %get facets from left looking at right (based on point only)
    overlapFacetsInds = find(intriangulation(fullPts, bFace, intersectionPoint(:,:,1)));

    %calculate combined FOV from indivdiual FOVs
    bincoularFOV = intersect(find(inFOVRight), find(inFOVLeft));
%     combinedFov = unique([find(inFOVLeft)' find(inFOVRight)']');
    outOfFov = intersect(find(inFOVLeft == 0), find(inFOVRight == 0));
    combinedFov = ceil((inFOVLeft + inFOVRight)/2);
    
    fprintf('monocular FOV: %.2f\n', sum(inFOVLeft)/size(headSphere,1)*100)
    fprintf('total FOV: %.2f\n', sum(combinedFov)/size(headSphere,1)*100)
    fprintf('binocular FOV: %.2f\n', length(bincoularFOV)/size(headSphere,1)*100)

    %calculate borders around FOV - use plotSphereR for later plotting
    [leftMapResults, numLeftReg] = calculateRegionsInMap( headSphere*plotSphereR, inFOVLeft, minDst*plotSphereR*4 );
    [rightMapResults, numRightReg] = calculateRegionsInMap( headSphere*plotSphereR, inFOVRight, minDst*plotSphereR*4 );
    inBino = zeros(length(inFOVLeft),1); inBino(bincoularFOV) = 1;
    [binoMapResults, numBinoReg] = calculateRegionsInMap( headSphere*plotSphereR, inBino, minDst*plotSphereR*4 );
       
    %- calculate vornoi diagram of eye and mapping to equirectangular image

    %firstly need coords for image
    imAz = (-180:360/(imageHeight*2-1):180)*-1;
    imEl = (-90:180/(imageHeight-1):90)*-1;

    imCoordsAz = zeros(imageHeight*2,imageHeight);
    imCoordsEl = zeros(imageHeight*2,imageHeight);
    for i = 1:imageHeight
        imCoordsAz(:,i) = imAz;
    end
    for i = 1:imageHeight*2
        imCoordsEl(i,:) = imEl;
    end

    imX = cos(imCoordsEl(:)/180*pi).*sin(imCoordsAz(:)/180*pi);
    imZ = -cos(imCoordsEl(:)/180*pi).*cos(imCoordsAz(:)/180*pi);
    imY = -sin(imCoordsEl(:)/180*pi);

    %calculate azimuth and elevation points on orignal head sphere
    
    %%%check az and el mapped correctly
    sphereEls = acos(headSphere(:,2))/pi*180-90;
    sphereAzs = atan2(headSphere(:,1),-headSphere(:,3))/pi*180;

    %get elvation range at zero degrees azimuth
    ZeroAzLineFullFrontInds = find(combinedFov & (sphereAzs == 0));
    ZeroAzLineLeftFrontInds = find(inFOVLeft & (sphereAzs == 0));
    ZeroAzLineFullRearTopInds = find(combinedFov & (abs(sphereAzs)== 180) & sphereEls > 0);
    ZeroAzLineLeftRearTopInds = find(inFOVLeft & (abs(sphereAzs)== 180) & sphereEls > 0);
    ZeroAzLineFullRearBotInds = find(combinedFov & (abs(sphereAzs)== 180) & sphereEls < 0);
    ZeroAzLineLeftRearBotInds = find(inFOVLeft & (abs(sphereAzs)== 180) & sphereEls < 0);
    
    if isempty(ZeroAzLineFullRearTopInds)
       maxElFull = max(sphereEls(ZeroAzLineFullFrontInds));
    else
       maxElFull = min(sphereEls(ZeroAzLineFullRearTopInds))+90;
    end
    if isempty(ZeroAzLineLeftRearTopInds)
       maxElLeft = max(sphereEls(ZeroAzLineLeftFrontInds));
    else
       maxElLeft = min(sphereEls(ZeroAzLineLeftRearTopInds))+90;
    end
    
    if isempty(ZeroAzLineFullRearBotInds)
       minElFull = min(sphereEls(ZeroAzLineFullFrontInds));
    else
       minElFull = max(sphereEls(ZeroAzLineFullRearBotInds))-90;
    end
    if isempty(ZeroAzLineLeftRearBotInds)
       minElLeft = min(sphereEls(ZeroAzLineLeftFrontInds));
    else
       minElLeft = max(sphereEls(ZeroAzLineLeftRearBotInds))-90;
    end
    
    ZeroElLineFullInds = find(combinedFov & (sphereEls == 0));
    ZeroElLineLeftInds = find(inFOVLeft & (sphereEls == 0));
    
    maxAzLeft = max(sphereAzs(ZeroElLineLeftInds));
    minAzLeft = min(sphereAzs(ZeroElLineLeftInds));
    maxAzFull = max(sphereAzs(ZeroElLineFullInds));
    minAzFull = min(sphereAzs(ZeroElLineFullInds));
    
    fprintf('Full: az range %.0f - %.0f (%.0f), el range %.0f - %.0f (%.0f)\n', minAzFull, maxAzFull, maxAzFull-minAzFull, minElFull, maxElFull, maxElFull-minElFull);
    fprintf('Left: az range %.0f - %.0f (%.0f), el range %.0f - %.0f (%.0f)\n', minAzLeft, maxAzLeft, maxAzLeft-minAzLeft, minElLeft, maxElLeft, maxElLeft-minElLeft);
    
    %%%if rear goes above zero degrees, determine manually from image using 
    %90-xxx/500*180 %for front
    %-(90*3-xxx/500*180) %for back
    
    %place vectors from original equirectangular image into voronoi mapping fit onto sphere
    facetAllocImageTemp = zeros(length(imAz), length(imEl)); 
    facetAllocImageTemp(:) = assignSubsToVoroniDiagram(intersectionPoint(:,:,1)',imX, imY, imZ );

    %interpolate which points on image are in left eye fov
    fovInterpolant = scatteredInterpolant( sphereAzs, sphereEls, inFOVLeft, 'nearest','nearest');
    fovImageLeft = zeros(length(imAz), length(imEl));
    fovImageLeft(:) = round(fovInterpolant(imCoordsAz(:), imCoordsEl(:)));

    %interpolate any bits of FOV that got left out of voroni label...
    misingINDs = find(fovImageLeft(:) == 1 & facetAllocImageTemp(:) == 0);
    goodINDs = find(fovImageLeft(:) == 1 & facetAllocImageTemp(:) > 1);
    fovInterpolant = scatteredInterpolant(imCoordsAz(goodINDs), imCoordsEl(goodINDs), facetAllocImageTemp(goodINDs), 'nearest', 'nearest');    
    facetAllocImageTemp(misingINDs) = round(fovInterpolant(imCoordsAz(misingINDs), imCoordsEl(misingINDs)));

%     figure; imshow_w_by_h(facetAllocImageTemp/50)
%     figure; imshow_w_by_h(fovImageLeft)
%     figure; imshow_w_by_h(imCoordsAz)
%     figure; imshow_w_by_h(imCoordsEl)
     
    %get coordinates for each facets, masked by field of view
    voroniCordsLeft = cell(totalPaths,1);
    for i = 1:length(unique(facetAllocImageTemp))
        voroniCordsLeft{i} = find(facetAllocImageTemp.*fovImageLeft == i);
    end

    %do same procedure for vectors from both
    facetAllocImageTemp = zeros(length(imAz), length(imEl)); 
    facetAllocImageTemp(:) = assignSubsToVoroniDiagram([intersectionPoint(:,:,1)' intersectionPoint(:,:,2)'],imX, imY, imZ );

    fovInterpolant = scatteredInterpolant( sphereAzs, sphereEls, combinedFov, 'nearest','nearest');
    fovImageBoth = zeros(length(imAz), length(imEl));
    fovImageBoth(:) = round(fovInterpolant(imCoordsAz(:), imCoordsEl(:)));

    %interpt missing
    misingINDs = find(fovImageBoth(:) == 1 & facetAllocImageTemp(:) == 0);
    goodINDs = find(fovImageBoth(:) == 1 & facetAllocImageTemp(:) > 1);
    fovInterpolant = scatteredInterpolant(imCoordsAz(goodINDs), imCoordsEl(goodINDs), facetAllocImageTemp(goodINDs), 'nearest', 'nearest');    
    facetAllocImageTemp(misingINDs) = round(fovInterpolant(imCoordsAz(misingINDs), imCoordsEl(misingINDs)));
    
    voroniCordsBoth = cell(totalPaths*2,1);
    for i = 1:length(unique(facetAllocImageTemp))
        voroniCordsBoth{i,1} = find(facetAllocImageTemp.*fovImageBoth == i);
    end
    
    %parse lines for 2D image...
    if useLargestRegion
        useBinoRegTo = 1;
        useRightRegTo = 1;
        binoImLines = cell(useBinoRegTo, 2);
        rightImLines = cell(useRightRegTo, 2);
        leftImLines = cell(useRightRegTo, 2);
    else 
        useBinoRegTo = numBinoReg;
        useRightRegTo = numRightReg;
        binoImLines = cell(useBinoRegTo, 2);
        rightImLines = cell(useRightRegTo, 2);
        leftImLines = cell(useRightRegTo, 2);
    end

    topInd = find(sphereEls == 90);
    bottomInd = find(sphereEls == -90);
    
    for i = 1:useBinoRegTo
        tempVerts = binoMapResults{i,1}/plotSphereR;
        lineEls = acos(tempVerts(:,2))/pi*180-90;
        lineAzs = atan2(tempVerts(:,1),-tempVerts(:,3))/pi*180;
    
        %assign to x and y coords matching el and az map
        xLeft = zeros(length(lineEls),1);
        yLeft = zeros(length(lineEls),1);
        for j = 1:length(lineEls)
            [~, xLeft(j)] = min(abs(imAz-lineAzs(j)));
            [~, yLeft(j)] = min(abs(imEl-lineEls(j)));
        end

          %solving wrapping problems
        if sum(binoMapResults{i,3} == topInd)
            tInd = find(abs(diff(xLeft)) > 50); %remove crappy jumps first
            xLeft(tInd) = [];
            yLeft(tInd) = [];
            %if so find point on left of linkage
            [~, leftLink] = min(xLeft);
            [~, rightLink] = max(xLeft);

            if rightLink > leftLink
                yLeft = [yLeft(1:leftLink)' yLeft(leftLink) 1 1 yLeft(rightLink) yLeft(rightLink:end)'];    
                xLeft = [xLeft(1:leftLink)' 1 1 imageHeight*2 imageHeight*2 xLeft(rightLink:end)'];
            else
                yLeft = [yLeft(1:rightLink)' yLeft(rightLink) 1 1 yLeft(leftLink) yLeft(leftLink:end)'];    
                xLeft = [xLeft(1:rightLink)' imageHeight*2 imageHeight*2 1 1 xLeft(leftLink:end)'];
            end
        end
        
        if sum(binoMapResults{i,3} == bottomInd)
            tInd = find(abs(diff(xLeft)) > 50); %remove crappy jumps first
            xLeft(tInd) = [];
            yLeft(tInd) = [];
            %find links
            [~, leftLink] = min(xLeft);
            [~, rightLink] = max(xLeft);
            if rightLink > leftLink
                yLeft = [yLeft(1:leftLink)' yLeft(leftLink) imageHeight imageHeight yLeft(rightLink) yLeft(rightLink:end)'];    
                xLeft = [xLeft(1:leftLink)' 1 1 imageHeight*2 imageHeight*2 xLeft(rightLink:end)'];
            else
                yLeft = [yLeft(1:rightLink)' yLeft(rightLink) imageHeight imageHeight yLeft(leftLink) yLeft(leftLink:end)'];    
                xLeft = [xLeft(1:rightLink)' imageHeight*2 imageHeight*2 1 1 xLeft(leftLink:end)'];
            end
        end
        
        binoImLines{i,1} = xLeft;
        binoImLines{i,2} = yLeft;
    end
    
    for i = 1:useRightRegTo
        tempVerts = rightMapResults{i,1}/plotSphereR;
        lineEls = acos(tempVerts(:,2))/pi*180-90;
        lineAzs = atan2(tempVerts(:,1),-tempVerts(:,3))/pi*180;
    
        %assign to x and y coords matching el and az map
        xLeft = zeros(length(lineEls),1);
        yLeft = zeros(length(lineEls),1);
        for j = 1:length(lineEls)
            [~, xLeft(j)] = min(abs(imAz-lineAzs(j)));
            [~, yLeft(j)] = min(abs(imEl-lineEls(j)));
        end
    
         %solving wrapping problems
        if sum(rightMapResults{i,3} == topInd)
            tInd = find(abs(diff(xLeft)) > 50); %remove crappy jumps first
            xLeft(tInd) = [];
            yLeft(tInd) = [];
            %if so find point on left of linkage
            [~, leftLink] = min(xLeft);
            [~, rightLink] = max(xLeft);

            if rightLink > leftLink
                yLeft = [yLeft(1:leftLink)' yLeft(leftLink) 1 1 yLeft(rightLink) yLeft(rightLink:end)'];    
                xLeft = [xLeft(1:leftLink)' 1 1 imageHeight*2 imageHeight*2 xLeft(rightLink:end)'];
            else
                yLeft = [yLeft(1:rightLink)' yLeft(rightLink) 1 1 yLeft(leftLink) yLeft(leftLink:end)'];    
                xLeft = [xLeft(1:rightLink)' imageHeight*2 imageHeight*2 1 1 xLeft(leftLink:end)'];
            end
        end
        
        if sum(rightMapResults{i,3} == bottomInd)
            tInd = find(abs(diff(xLeft)) > 50); %remove crappy jumps first
            xLeft(tInd) = [];
            yLeft(tInd) = [];
            %find links
            [~, leftLink] = min(xLeft);
            [~, rightLink] = max(xLeft);
            if rightLink > leftLink
                yLeft = [yLeft(1:leftLink)' yLeft(leftLink) imageHeight imageHeight yLeft(rightLink) yLeft(rightLink:end)'];    
                xLeft = [xLeft(1:leftLink)' 1 1 imageHeight*2 imageHeight*2 xLeft(rightLink:end)'];
            else
                yLeft = [yLeft(1:rightLink)' yLeft(rightLink) imageHeight imageHeight yLeft(leftLink) yLeft(leftLink:end)'];    
                xLeft = [xLeft(1:rightLink)' imageHeight*2 imageHeight*2 1 1 xLeft(leftLink:end)'];
            end
        end
        rightImLines{i,1} = xLeft;
        rightImLines{i,2} = yLeft;
        
        %do for left as well
        tempVerts = leftMapResults{i,1}/plotSphereR;
        lineEls = acos(tempVerts(:,2))/pi*180-90;
        lineAzs = atan2(tempVerts(:,1),-tempVerts(:,3))/pi*180;
    
        %assign to x and y coords matching el and az map
        xLeft = zeros(length(lineEls),1);
        yLeft = zeros(length(lineEls),1);
        for j = 1:length(lineEls)
            [~, xLeft(j)] = min(abs(imAz-lineAzs(j)));
            [~, yLeft(j)] = min(abs(imEl-lineEls(j)));
        end
    
         %solving wrapping problems
        if sum(leftMapResults{i,3} == topInd)
            tInd = find(abs(diff(xLeft)) > 50); %remove crappy jumps first
            xLeft(tInd) = [];
            yLeft(tInd) = [];
            %if so find point on left of linkage
            [~, leftLink] = min(xLeft);
            [~, rightLink] = max(xLeft);

            if rightLink > leftLink
                yLeft = [yLeft(1:leftLink)' yLeft(leftLink) 1 1 yLeft(rightLink) yLeft(rightLink:end)'];    
                xLeft = [xLeft(1:leftLink)' 1 1 imageHeight*2 imageHeight*2 xLeft(rightLink:end)'];
            else
                yLeft = [yLeft(1:rightLink)' yLeft(rightLink) 1 1 yLeft(leftLink) yLeft(leftLink:end)'];    
                xLeft = [xLeft(1:rightLink)' imageHeight*2 imageHeight*2 1 1 xLeft(leftLink:end)'];
            end
        end
        
        if sum(rightMapResults{i,3} == bottomInd)
            tInd = find(abs(diff(xLeft)) > 50); %remove crappy jumps first
            xLeft(tInd) = [];
            yLeft(tInd) = [];
            %find links
            [~, leftLink] = min(xLeft);
            [~, rightLink] = max(xLeft);
            if rightLink > leftLink
                yLeft = [yLeft(1:leftLink)' yLeft(leftLink) imageHeight imageHeight yLeft(rightLink) yLeft(rightLink:end)'];    
                xLeft = [xLeft(1:leftLink)' 1 1 imageHeight*2 imageHeight*2 xLeft(rightLink:end)'];
            else
                yLeft = [yLeft(1:rightLink)' yLeft(rightLink) imageHeight imageHeight yLeft(leftLink) yLeft(leftLink:end)'];    
                xLeft = [xLeft(1:rightLink)' imageHeight*2 imageHeight*2 1 1 xLeft(leftLink:end)'];
            end
        end
        leftImLines{i,1} = xLeft;
        leftImLines{i,2} = yLeft;
    end
end

%--- plotting
if doPlotting
    %-plot fields of view on sphere
    r = 0:0.01:2*pi;
    frontalRing = [sin(r)'*plotSphereR, cos(r)'*plotSphereR, zeros(length(r),1)];
    sagitalRing = [zeros(length(r),1), sin(r)'*plotSphereR, cos(r)'*plotSphereR];
    axialRing = [sin(r)'*plotSphereR, zeros(length(r),1), cos(r)'*plotSphereR];
    ring45 = [sin(r)'*plotSphereR*0.7071, ones(length(r),1)*-0.7071*plotSphereR, cos(r)'*plotSphereR*0.7071];
    ring45N = [sin(r)'*plotSphereR*0.7071, ones(length(r),1)*0.7071*plotSphereR, cos(r)'*plotSphereR*0.7071];
   
    if bgColor
        bgCol = [1 1 1];
        textCol = [0 0 0];
    else
        bgCol = [0 0 0];
        textCol = [1 1 1];
    end
    
    vFOnSphereF = figure; hold on
    patch(struct('faces', F,'vertices',Vshifted),'FaceColor', [0.2 0.8 0.2], 'EdgeColor', 'none','FaceNormals',eyeNormsShiftedRot,...
        'SpecularStrength', 0.4, 'DiffuseStrength', 0.7, 'AmbientStrength', 0.7, 'SpecularExponent', 2); %[0 0 0]
    if ~isempty(headSTLFile)
        patch(struct('faces', FHead,'vertices',VHeadShifted),'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none','FaceNormals',HeadNormsRot,...
             'SpecularStrength', 0.3, 'DiffuseStrength', 0.7, 'AmbientStrength', 0.7, 'SpecularExponent', 2); %[0 0 0]); %[0 0 0]
    end
     axis equal; axis vis3d; axis off
    set(gcf, 'Color',bgCol); set(gca,'clipping','off');

    sphereCols = ones(size(TR.X,1),3); 
    sphereCols(find(inFOVLeft),1) = 0; sphereCols(find(inFOVLeft),2) = 1; sphereCols(find(inFOVLeft),3) = 0;
    sphereCols(find(inFOVRight),1) = 0; sphereCols(find(inFOVRight),2) = 0; sphereCols(find(inFOVRight),3) = 1;
    sphereCols(bincoularFOV,1) = 0; sphereCols(bincoularFOV,2) = 1; sphereCols(bincoularFOV,3) = 1;
    patch(struct('faces', TR.Triangulation,'vertices',TR.X*plotSphereR),'FaceVertexCData', sphereCols,'EdgeColor', 'none', 'EdgeAlpha', 0,...
        'FaceAlpha', 0.6, 'FaceColor', 'interp', 'SpecularStrength', 0, 'DiffuseStrength', 0.8, 'AmbientStrength', 0.8); 
    
    scatter3(headSphere(closestPoint,1)*plotSphereR, headSphere(closestPoint,2)*plotSphereR, headSphere(closestPoint,3)*plotSphereR, 10,...
       'MarkerEdgeColor','k', 'MarkerFaceColor',[0 1 0])

   scatter3(headSphere(closestPointRight,1)*plotSphereR, headSphere(closestPointRight,2)*plotSphereR, headSphere(closestPointRight,3)*plotSphereR, 10,...
       'MarkerEdgeColor','k', 'MarkerFaceColor',[0 0 1])
   
    plot3(frontalRing(:,1), frontalRing(:,2), frontalRing(:,3), '-', 'linewidth', 3, 'color', [0.7 0.7 0.7]);  
    plot3(sagitalRing(:,1), sagitalRing(:,2), sagitalRing(:,3), '-', 'linewidth', 3, 'color', [0.7 0.7 0.7]);  
    plot3(axialRing(:,1), axialRing(:,2), axialRing(:,3), '-', 'linewidth', 3, 'color', [0.7 0.7 0.7]);  
    plot3(ring45(:,1), ring45(:,2), ring45(:,3), '-', 'linewidth', 3, 'color', [0.7 0.7 0.7]);  
    plot3(ring45N(:,1), ring45N(:,2), ring45N(:,3), '-', 'linewidth', 3, 'color', [0.7 0.7 0.7]); 
    
    hM1 = mArrow3([0 0  0], [ plotSphereR*1.4 0 0], 'color', 'g', 'FaceLighting','none');
    hM2 = mArrow3([0 0  0], [ 0 plotSphereR*1.4*-1 0], 'color', 'r', 'FaceLighting','none');
    hM3 = mArrow3([0 0  0], [ 0 0 plotSphereR*1.4*-1], 'color', 'c', 'FaceLighting','none');
    if plotRightEyeLine
        for i = 1: numLeftReg
            tempVerts = leftMapResults{i,1};
            plot3(tempVerts(:,1), tempVerts(:,2), tempVerts(:,3), '-', 'color', 'g', 'linewidth', 2);
        end

        for i = 1:useRightRegTo
            tempVerts = rightMapResults{i,1};
            plot3(tempVerts(:,1), tempVerts(:,2), tempVerts(:,3), '-', 'color', 'b', 'linewidth', 2);
        end
    end
    
    if plotBinoLine
        for i = 1:useBinoRegTo
            tempVerts = binoMapResults{i,1};
            plot3(tempVerts(:,1), tempVerts(:,2), tempVerts(:,3), '-', 'color', 'c', 'linewidth', 2);
        end
    end
    
    if displayTitles
       title(sprintf('Visual fields of both eyes projected onto a sphere\n(cyan indicates binocular overlap in all plots)'), 'color', textCol); 
    end
    if autoAlignAndIllum
        view(AlignAngle);
        camlight('headlight'); lighting gouraud;
    end
    set(gcf, 'Position', [50 50 1000 1000]);
           
    %- plot things on eye and world for each parameter
    %- lens d
    plotLensDEyeF = figure;
    patch(struct('faces', F,'vertices',Vshifted),'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none','FaceNormals',eyeNormsShiftedRot,...
        'SpecularStrength', 0.4, 'DiffuseStrength', 0.7, 'AmbientStrength', 0.7, 'SpecularExponent', 2); %[0 0 0]
    if ~isempty(headSTLFile)
        patch(struct('faces', FHead,'vertices',VHeadShifted),'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none','FaceNormals',HeadNormsRot,...
             'SpecularStrength', 0.3, 'DiffuseStrength', 0.7, 'AmbientStrength', 0.7, 'SpecularExponent', 2); %[0 0 0]); %[0 0 0]
    end
     axis equal; axis vis3d; axis off; hold on; set(gca,'clipping','off');
     
    if exist('xDLower') & exist('xIOLower') & isempty(numColBins)
        numColBins = length(xIOLower);
        usePreloaded = 1;
    else
        if isempty(numColBins)
            numColBins = 10;
        end
        usePreloaded = 0;
    end
    
    colsLD = (copper(numColBins+1));
    colsLD(1,:) = []; %remove 1st as black...
    colsLD = [colsLD(:,2), colsLD(:,3), colsLD(:,1)];
    colIds = zeros(length(lensDiamteres),1);
    
    colImg = zeros(imageHeight*2, imageHeight,3);
    lensDUWImage = zeros(imageHeight*2, imageHeight);
    if sum(bgColor)
        colImg(:) = 1;
    end
    if exist('xDLower') & usePreloaded
        for i = 1:numColBins-1
            if i == 1
                inds = find(lensDiamteres >= 0 & lensDiamteres < xDLower(i+1));
            elseif i == numColBins-1
                inds = find(lensDiamteres >= xDLower(i) & lensDiamteres < 1000);
            else
                inds = find(lensDiamteres >= xDLower(i) & lensDiamteres < xDLower(i+1));
            end
            colIds(inds) = i;
            for j = 1:length(inds)
                [sX, sY] = ind2sub([imageHeight*2, imageHeight],voroniCordsLeft{inds(j)});
                lensDUWImage(voroniCordsLeft{inds(j)}) = mean(xDLower([i i+1]));
                for k = 1:length(sX)
                    colImg(sX(k), sY(k),:) = colsLD(i,:);
                end
            end
        end
    else
       %break histogram up into equal distribution ranges
       [n, x] = hist(lensDiamteres,100);
       cumSAngle = cumsum(n)/sum(n);
       xDiaLowerTemp = zeros(numColBins+1,1);
       xDiaLowerTemp(1) = x(1)-1; xDiaLowerTemp(end) = x(end)+1;
       for i = 1:numColBins-1
             temp = find(cumSAngle > i/numColBins);
             xDiaLowerTemp(i+1) = x(temp(1));
       end
       
       for i = 1:numColBins
            inds = find(lensDiamteres >= xDiaLowerTemp(i) & lensDiamteres < xDiaLowerTemp(i+1));
            colIds(inds) = i;
            for j = 1:length(inds)
                [sX, sY] = ind2sub([imageHeight*2, imageHeight],voroniCordsLeft{inds(j)});
                lensDUWImage(voroniCordsLeft{inds(j)}) = mean(xDiaLowerTemp([i i+1]));
                for k = 1:length(sX)
                    colImg(sX(k), sY(k),:) = colsLD(i,:);
                end
            end
       end
    end
    
    for i = 1:totalPaths
        tempVerts = eachPatch{i,1};
        plot3(tempVerts(:,1), tempVerts(:,2), tempVerts(:,3), '.', 'color', colsLD(colIds(i),:));
        
        if sum(overlapFacetsInds == i) & ~isempty(mSize)
            plot3(patchTips(i,1,1), patchTips(i,2,1), patchTips(i,3,1), 'c.', 'markersize',mSize);
        end
    end
    set(gcf, 'Color',bgCol); 
    
    if displayTitles
       title(sprintf('Lens diameters (D) displayed on eye\n(Color bar used for associated unwrapped image)'), 'color', textCol); 
    end
    
    if autoAlignAndIllum
        view(AlignAngle);
        camlight('headlight'); lighting gouraud;
    end
   set(gcf, 'Position', [50 50 800 800]);
           
    DColorbar = figure;
    labels = cell(numColBins,1);
    if exist('xDLower') & usePreloaded
        labs2use = xDLower;
    else
        labs2use = xDiaLowerTemp;
    end
    %make into linear representation of color values
%     tempColVals = zeros(101,3);
%     tempColVals(1:10,1) = colsLD(1,1); tempColVals(1:10,2) = colsLD(1,2); tempColVals(1:10,3) = colsLD(1,3);
%     tempColVals(91:101,1) = colsLD(end,1); tempColVals(91:101,2) = colsLD(end,2); tempColVals(91:101,3) = colsLD(end,3);
%     
%     cbBarSteps = zeros(numColBins,1);
%     cbBarSteps(1) = 5;
%     cbBarSteps(end) = 96;
%     minVal = mean(labs2use(1:2)); maxVal =  mean(labs2use(end-1:end));
%     
%     for i = 2:numColBins-1
%         centStep = mean(labs2use(i:i+1));
%         cbBarSteps(i) = (centStep-minVal)/(maxVal-minVal)*(cbBarSteps(end)-cbBarSteps(1))+cbBarSteps(1);
%     end
    colormap(colsLD); hcb = colorbar;
    
    for i = 1:numColBins
        labels{i} = sprintf('>%.1f', labs2use(i));
%         if i == 1
%             labels{i} = sprintf('<%.1f', labs2use(i+1));
%         else if i == numColBins
%                 labels{i} = sprintf('>%.1f', labs2use(i));
%             else
%                 labels{i} = sprintf('%.1f - %.1f', labs2use(i), labs2use(i+1));
%         end; end
    end
    set(hcb,'Color', textCol, 'TickDirection','out',...
        'Ticks', (0:1/(numColBins):1-1/(numColBins)) + (1/(numColBins))/2,...
        'TickLabels', labels);
    ylabel(hcb,sprintf('Lens Diameter (%sm)', '\mu'));
    set(gcf, 'Color',bgCol); 
    
    %do world image for diameters
    diamUWFig = figure; 
    imshow_w_by_h(colImg); hold on
    
    if plotBinoLine
        for i = 1:useBinoRegTo
            plot(binoImLines{i,1}, binoImLines{i,2}, 'c-', 'linewidth',2);
        end
    end
    
    if plotRightEyeLine
        for i = 1:useRightRegTo
            plot(rightImLines{i,1}, rightImLines{i,2}, 'b-', 'linewidth',2);
        end
    end
    
    if displayTitles
       title(sprintf('Lens diameters (D) displayed on unwrapped image')); 
    end
    if displayAxisLabelsOnUnWrappedImages
       xlabel('Azimuth(^o)');
       ylabel('Elevation(^o)');
       axis on
       set(gca, 'Xtick', [1 imageHeight/2 imageHeight imageHeight*3/2 imageHeight*2], 'XTickLabel', {'\pm180', '-90', '0', '90', '\pm180'},...
           'Ytick', [1 imageHeight*0.25 imageHeight/2 imageHeight*0.75 imageHeight], 'YTickLabel', {'90', '45', '0', '-45', '-90'},'box','off');
    end
    
    %- interommatidial angles
    plotIOEyeF = figure;
    patch(struct('faces', F,'vertices',Vshifted),'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none','FaceNormals',eyeNormsShiftedRot,...
        'SpecularStrength', 0.4, 'DiffuseStrength', 0.7, 'AmbientStrength', 0.7, 'SpecularExponent', 2); %[0 0 0]
    
    if ~isempty(headSTLFile)
        patch(struct('faces', FHead,'vertices',VHeadShifted),'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none','FaceNormals',HeadNormsRot,...
             'SpecularStrength', 0.3, 'DiffuseStrength', 0.7, 'AmbientStrength', 0.7, 'SpecularExponent', 2); %[0 0 0]); %[0 0 0]
    end
    axis equal; axis vis3d; axis off; set(gca,'clipping','off'); hold on
     
    colsIO = flipud(copper(numColBins+1));
    colsIO(end,:) = []; %remove last as black
    colIds = zeros(length(interOAngle),1);
    
    colImg = zeros(imageHeight*2, imageHeight,3);
    interOUWImage = zeros(imageHeight*2, imageHeight);
    if sum(bgColor)
        colImg(:) = 1;
    end
    if exist('xIOLower') & usePreloaded
        for i = 1:numColBins-1
            if i == 1
                inds = find(interOAngle >= 0 & interOAngle < xIOLower(i+1));
            elseif i == numColBins-1
                inds = find(interOAngle >= xIOLower(i) & interOAngle < 1000);
            else
                inds = find(interOAngle >= xIOLower(i) & interOAngle < xIOLower(i+1));
            end
            colIds(inds) = i;
            for j = 1:length(inds)
                [sX, sY] = ind2sub([imageHeight*2, imageHeight],voroniCordsLeft{inds(j)});
                interOUWImage(voroniCordsLeft{inds(j)}) = mean(xIOLower([i i+1]));
                for k = 1:length(sX)
                    colImg(sX(k), sY(k),:) = colsIO(i,:);
                end
            end
        end
    else
       %break histogram up into equal distribution ranges
       [n, x] = hist(interOAngle,100);
       cumSAngle = cumsum(n)/sum(n);
       xPhiLowerTemp = zeros(numColBins+1,1);
       xPhiLowerTemp(1) = x(1)-1; xPhiLowerTemp(end) = x(end)+1;
       for i = 1:numColBins-1
             temp = find(cumSAngle > i/numColBins);
             xPhiLowerTemp(i+1) = x(temp(1));
       end
       
       for i = 1:numColBins
            inds = find(interOAngle >= xPhiLowerTemp(i) & interOAngle < xPhiLowerTemp(i+1));
            colIds(inds) = i;
            for j = 1:length(inds)
                [sX, sY] = ind2sub([imageHeight*2, imageHeight],voroniCordsLeft{inds(j)});
                interOUWImage(voroniCordsLeft{inds(j)}) = mean(xPhiLowerTemp([i i+1]));
                for k = 1:length(sX)
                    colImg(sX(k), sY(k),:) = colsIO(i,:);
                end
            end
       end
    end
    
    for i = 1:totalPaths
        tempVerts = eachPatch{i,1};
        plot3(tempVerts(:,1), tempVerts(:,2), tempVerts(:,3), '.', 'color', colsIO(colIds(i),:));
        
        if sum(overlapFacetsInds == i) & ~isempty(mSize)
            plot3(patchTips(i,1,1), patchTips(i,2,1), patchTips(i,3,1), 'c.', 'markersize',mSize);
        end
    end

    set(gcf, 'Color',bgCol); 
    
    if displayTitles
       title(sprintf('Interommatidial angles (%s) displayed on eye\n(Color bar used for associated unwrapped image)', '\Delta\phi'), 'color', textCol); 
    end
    
    if autoAlignAndIllum
        view(AlignAngle);
        camlight('headlight'); lighting gouraud;
    end
    set(gcf, 'Position', [50 50 800 800]);
    
    IOColourbarFig = figure;
    colormap(colsIO); hcb = colorbar;
    labels = cell(numColBins,1);
    if exist('xIOLower') & usePreloaded
        labs2use = xIOLower;
    else
        labs2use = xPhiLowerTemp;
    end
    for i = 1:numColBins
        labels{i} = sprintf('>%.1f', labs2use(i));
%         if i == 1
%             labels{i} = sprintf('<%.1f', labs2use(i+1));
%         else if i == numColBins
%                 labels{i} = sprintf('>%.1f', labs2use(i));
%             else
%                 labels{i} = sprintf('%.1f - %.1f', labs2use(i), labs2use(i+1));
%         end; end
    end
    set(hcb,'Color', textCol, 'TickDirection','out',...
        'Ticks', (0:1/(numColBins):1-1/(numColBins)) + (1/(numColBins))/2,...
        'TickLabels', labels);
    ylabel(hcb,'Interommatidial Angle (^o)');
    set(gcf, 'Color',bgCol); 
    
    %do world image for angles
    interOUWFig = figure; 
    imshow_w_by_h(colImg); hold on
    
    if plotBinoLine
        for i = 1:useBinoRegTo
            plot(binoImLines{i,1}, binoImLines{i,2}, 'c-', 'linewidth',2);
        end
    end
    
    if plotRightEyeLine
        for i = 1:useRightRegTo
            plot(rightImLines{i,1}, rightImLines{i,2}, 'b-', 'linewidth',2);
        end
    end
    
    if displayTitles
       title(sprintf('Interommatidial angles (%s) displayed on unwrapped image', '\Delta\phi')); 
    end
    if displayAxisLabelsOnUnWrappedImages
       xlabel('Azimuth(^o)');
       ylabel('Elevation(^o)');
       axis on
       set(gca, 'Xtick', [1 imageHeight/2 imageHeight imageHeight*3/2 imageHeight*2], 'XTickLabel', {'\pm180', '-90', '0', '90', '\pm180'},...
           'Ytick', [1 imageHeight*0.25 imageHeight/2 imageHeight*0.75 imageHeight], 'YTickLabel', {'90', '45', '0', '-45', '-90'},'box','off');
    end
    
    %-eye parameter
    plotPEyeF = figure;
    patch(struct('faces', F,'vertices',Vshifted),'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none','FaceNormals',eyeNormsShiftedRot,...
        'SpecularStrength', 0.4, 'DiffuseStrength', 0.7, 'AmbientStrength', 0.7, 'SpecularExponent', 2); %[0 0 0]
    
    if ~isempty(headSTLFile)
        patch(struct('faces', FHead,'vertices',VHeadShifted),'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none','FaceNormals',HeadNormsRot,...
             'SpecularStrength', 0.3, 'DiffuseStrength', 0.7, 'AmbientStrength', 0.7, 'SpecularExponent', 2); %[0 0 0]); %[0 0 0]
    end
    axis equal; axis vis3d; axis off; set(gca,'clipping','off'); hold on
     
    colsLP = (hot(numColBins+2));
    colsLP(end,:) = []; %remove end as too white...
    colsLP(end-1,:) = []; %remove 2nd last as too yellow...
    colIds = zeros(length(interOAngle),1);
    
    colImg = zeros(imageHeight*2, imageHeight,3);
    lensPUWImage = zeros(imageHeight*2, imageHeight);
    if sum(bgColor)
        colImg(:) = 1;
    end
    if exist('xPLower') & usePreloaded
        for i = 1:numColBins-1
            if i == 1
                inds = find(eyeParam >= 0 & eyeParam < xPLower(i+1));
            elseif i == numColBins-1
                inds = find(eyeParam >= xPLower(i) & eyeParam < 1000);
            else
                inds = find(eyeParam >= xPLower(i) & eyeParam < xPLower(i+1));
            end
            colIds(inds) = i;
            for j = 1:length(inds)
                [sX, sY] = ind2sub([imageHeight*2, imageHeight],voroniCordsLeft{inds(j)});
                lensPUWImage(voroniCordsLeft{inds(j)}) = mean(xPLower([i i+1]));
                for k = 1:length(sX)
                    colImg(sX(k), sY(k),:) = colsLP(i,:);
                end
            end
        end
    else
       %break histogram up into equal distribution ranges
       [n, x] = hist(eyeParam,100);
       cumSAngle = cumsum(n)/sum(n);
       xPLowerTemp = zeros(numColBins+1,1);
       xPLowerTemp(1) = x(1)-1; xPLowerTemp(end) = x(end)+1;
       for i = 1:numColBins-1
             temp = find(cumSAngle > i/numColBins);
             xPLowerTemp(i+1) = x(temp(1));
       end
       
       for i = 1:numColBins
            inds = find(eyeParam >= xPLowerTemp(i) & eyeParam < xPLowerTemp(i+1));
            colIds(inds) = i;
            for j = 1:length(inds)
                [sX, sY] = ind2sub([imageHeight*2, imageHeight],voroniCordsLeft{inds(j)});
                lensPUWImage(voroniCordsLeft{inds(j)}) = mean(xPLowerTemp([i i+1]));
                for k = 1:length(sX)
                    colImg(sX(k), sY(k),:) = colsLP(i,:);
                end
            end
       end
    end
    
    for i = 1:totalPaths
        tempVerts = eachPatch{i,1};
        plot3(tempVerts(:,1), tempVerts(:,2), tempVerts(:,3), '.', 'color', colsLP(colIds(i),:));
        
        if sum(overlapFacetsInds == i) & ~isempty(mSize)
            plot3(patchTips(i,1,1), patchTips(i,2,1), patchTips(i,3,1), 'c.', 'markersize',mSize);
        end
    end
  
    set(gcf, 'Color',bgCol); 
    
    if displayTitles
       title(sprintf('Eye parameter (P) displayed on eye\n(Color bar used for associated unwrapped image)'), 'color', textCol); 
    end
    
    if autoAlignAndIllum
        view(AlignAngle);
        camlight('headlight'); lighting gouraud;
    end
    set(gcf, 'Position', [50 50 800 800]);
           
    PColourbarFig = figure;
    colormap(colsLP); hcb = colorbar;
    labels = cell(numColBins,1);
    if exist('xPLower') & usePreloaded
        labs2use = xPLower;
    else
        labs2use = xPLowerTemp;
    end
    for i = 1:numColBins
        labels{i} = sprintf('>%.1f', labs2use(i));
%         if i == 1
%             labels{i} = sprintf('<%.1f', labs2use(i+1));
%         else if i == numColBins
%                 labels{i} = sprintf('>%.1f', labs2use(i));
%             else
%                 labels{i} = sprintf('%.1f - %.1f', labs2use(i), labs2use(i+1));
%         end; end
    end
    set(hcb,'Color', textCol, 'TickDirection','out',...
        'Ticks', (0:1/(numColBins):1-1/(numColBins)) + (1/(numColBins))/2,...
        'TickLabels', labels);
    ylabel(hcb,sprintf('Eye Parameter (%s.rad)', '\mu'));
    set(gcf, 'Color',bgCol); 
    
    %do world image for parameters
    lensPUWFig = figure; 
    imshow_w_by_h(colImg); hold on
      
    if plotBinoLine
        for i = 1:useBinoRegTo
            plot(binoImLines{i,1}, binoImLines{i,2}, 'c-', 'linewidth',2);
        end
    end
    
    if plotRightEyeLine
        for i = 1:useRightRegTo
            plot(rightImLines{i,1}, rightImLines{i,2}, 'b-', 'linewidth',2);
        end
    end

    if displayTitles
       title(sprintf('Eye parameter (P) displayed on unwrapped image')); 
    end
    if displayAxisLabelsOnUnWrappedImages
       xlabel('Azimuth(^o)');
       ylabel('Elevation(^o)');
       axis on
       set(gca, 'Xtick', [1 imageHeight/2 imageHeight imageHeight*3/2 imageHeight*2], 'XTickLabel', {'\pm180', '-90', '0', '90', '\pm180'},...
           'Ytick', [1 imageHeight*0.25 imageHeight/2 imageHeight*0.75 imageHeight], 'YTickLabel', {'90', '45', '0', '-45', '-90'},'box','off');
    end
    
    %- plot contour plots if desired
    lineWidths = fliplr(1:2/(numColBins-1):3); %helps with make sequantially drawn lines visible
    
    if ~isempty(dispContourOnDiameterImage)
       figure(diamUWFig);
       if ~isempty(find(dispContourOnDiameterImage == 2))
           vals = unique(interOUWImage(:));
           vals(1) = []; %should be equal length to numColBins
           zInds = find(interOUWImage == 0);
           interOUWImage(zInds) = vals(end); %do this to avoid zero border error
           for i = 1:length(vals)
             if i == length(vals)
                 interOUWImage(zInds) = 0;
             end  
             contour(interOUWImage', [vals(i)*0.9 vals(i)*1.1], 'ShowText','off', 'linewidth', lineWidths(i), 'color', colsIO(i,:));
           end
       end
       
       if ~isempty(find(dispContourOnDiameterImage == 3))
           vals = unique(lensPUWImage(:));
           vals(1) = []; %should be equal length to numColBins
           zInds = find(lensPUWImage == 0);
           lensPUWImage(zInds) = vals(end); %do this to avoid zero border error
           for i = 1:length(vals)
             if i == length(vals)
                 lensPUWImage(zInds) = 0;
             end  
             contour(lensPUWImage', [vals(i)*0.9 vals(i)*1.1], 'ShowText','off', 'linewidth', lineWidths(i), 'color', colsLP(i,:));
           end
       end
    end
    
     if ~isempty(dispContourOnInterOImage)
       figure(interOUWFig);
       if ~isempty(find(dispContourOnInterOImage == 1))
           vals = unique(lensDUWImage(:));
           vals(1) = []; %should be equal length to numColBins
           zInds = find(lensDUWImage == 0);
           lensDUWImage(zInds) = vals(end); %do this to avoid zero border error
           for i = 1:length(vals)
             if i == length(vals)
                 lensDUWImage(zInds) = 0;
             end  
             contour(lensDUWImage', [vals(i)*0.9 vals(i)*1.1], 'ShowText','off', 'linewidth', lineWidths(i), 'color', colsLD(i,:));
           end
       end
       
       if ~isempty(find(dispContourOnInterOImage == 3))
           vals = unique(lensPUWImage(:));
           vals(1) = []; %should be equal length to numColBins
           zInds = find(lensPUWImage == 0);
           lensPUWImage(zInds) = vals(end); %do this to avoid zero border error
           for i = 1:length(vals)
             if i == length(vals)
                 lensPUWImage(zInds) = 0;
             end  
             contour(lensPUWImage', [vals(i)*0.9 vals(i)*1.1], 'ShowText','off', 'linewidth', lineWidths(i), 'color', colsLP(i,:));
           end
       end
     end
    
     if ~isempty(dispContourOnLensPImage)
       figure(lensPUWFig);
       if ~isempty(find(dispContourOnLensPImage == 1))
           vals = unique(lensDUWImage(:));
           vals(1) = []; %should be equal length to numColBins
           zInds = find(lensDUWImage == 0);
           lensDUWImage(zInds) = vals(end); %do this to avoid zero border error
           for i = 1:length(vals)
             if i == length(vals)
                 lensDUWImage(zInds) = 0;
             end  
             contour(lensDUWImage', [vals(i)*0.9 vals(i)*1.1], 'ShowText','off', 'linewidth', lineWidths(i), 'color', colsLD(i,:));
           end
       end
       
       if ~isempty(find(dispContourOnLensPImage == 2))
           vals = unique(interOUWImage(:));
           vals(1) = []; %should be equal length to numColBins
           zInds = find(interOUWImage == 0);
           interOUWImage(zInds) = vals(end); %do this to avoid zero border error
           for i = 1:length(vals)
             if i == length(vals)
                 interOUWImage(zInds) = 0;
             end  
             contour(interOUWImage', [vals(i)*0.9 vals(i)*1.1], 'ShowText','off', 'linewidth', lineWidths(i), 'color', colsIO(i,:));
           end
       end
    end
    
    %- create images for simulation
    if simUserImage
        %create simulated image
        workingIm = imread_w_by_h(imageToGet,imType);
        workingIm = rgb2gray(workingIm);
        workingIm = double(imresize(workingIm,[imageHeight*2,imageHeight]))/255;
        if imgRotate ~= 0
           imgShift = (imgRotate/180*500);
           if imgShift < 1
                workingIm = [workingIm(end+imgShift:end,:)' (workingIm(1:end+imgShift,:))']'; 
           else
                workingIm = [workingIm(imgShift:end,:)' (workingIm(1:imgShift,:))']'; 
           end
        end
        if contrastAdj ~= 0
           workingIm = imadjust(workingIm,[contrastAdj 1-contrastAdj]); 
        end
%         workingIm = im2bw(workingIm,0.5);
%         figure;imshow_w_by_h(workingIm);
    else
        %create check image
        workingIm = zeros(length(imAz), length(imEl));
        azSteps = numChequers;
        elSteps = numChequers/2; %el should be half az for isotropic grid
        azBins = -180:360/azSteps:180;
        elBins = -90:180/elSteps:90;

        startSet = 0;
        for i = 1:length(azBins)-1
            lineSet = startSet;
           for j = 1:length(elBins)-1

               inds = find(imCoordsAz >= azBins(i) & imCoordsAz <= azBins(i+1) &...
                   imCoordsEl >= elBins(j) & imCoordsEl <= elBins(j+1));

               workingIm(inds) = lineSet;
               lineSet = ~lineSet;
           end
           startSet = ~startSet;
        end
    end
    
    %plot it with lines...
    linesOnOrigF = figure; imshow_w_by_h(workingIm); hold on
    if plotBinoLine
        for i = 1:useBinoRegTo
            plot(binoImLines{i,1}, binoImLines{i,2}, 'c-', 'linewidth',2);
        end
    end
    
    if plotRightEyeLine
        for i = 1: useRightRegTo
            plot(leftImLines{i,1}, leftImLines{i,2}, 'g-', 'linewidth',2);
        end

        for i = 1:useRightRegTo
            plot(rightImLines{i,1}, rightImLines{i,2}, 'b-', 'linewidth',2);
        end
    end
    
    if displayTitles
        if ~simUserImage
            title(sprintf('original unwrapped image for simulation %.0f^o checks with %.2f cyc/^o period', 360/numChequers, 1/(360/numChequers*2))); 
        else
            title('Original unwrapped image for simulation'); 
        end
    end
    if displayAxisLabelsOnUnWrappedImages
       xlabel('Azimuth(^o)');
       ylabel('Elevation(^o)');
       axis on
       set(gca, 'Xtick', [1 imageHeight/2 imageHeight imageHeight*3/2 imageHeight*2], 'XTickLabel', {'\pm180', '-90', '0', '90', '\pm180'},...
           'Ytick', [1 imageHeight*0.25 imageHeight/2 imageHeight*0.75 imageHeight], 'YTickLabel', {'-90', '-45', '0', '45', '90'},'box','off');
    end
    
    %- simulate gnatt view
    if sensitvityScaleByColSat
        satImg = zeros(length(imAz), length(imEl));     
    end
    newImage = zeros(length(imAz), length(imEl)); 
    if simBothEyes
       numToGet = size(voroniCordsBoth,1);
    else
       numToGet = size(voroniCordsLeft,1);
    end
    for i = 1:numToGet
        if i <= totalPaths
            ref = i;
            setC = 1;
        else
            ref = i-totalPaths;
            setC = 2;
        end
%         crdDist = 2*sin(acceptanceAngs(ref)/180*pi/2/2*3); %find distance of a cord given the acceptance angle
%         indsIncluded = find(sqrt((intersectionPoint(ref,1,setC)-imX).^2 + (intersectionPoint(ref,2,setC)-imY).^2 + (intersectionPoint(ref,3,setC)-imZ).^2) <= crdDist);

         deltaAng = 2*asin(sqrt((intersectionPoint(ref,1,setC)-imX).^2+(intersectionPoint(ref,2,setC)-imY).^2+...
                      (intersectionPoint(ref,3,setC)-imZ).^2)/2); 
         indsIncluded = find(deltaAng < acceptanceAngs(ref)/180*pi);
         
         anglesToOurPoint = deltaAng(indsIncluded);
         
        %anglesToOurPoint = zeros(length(indsIncluded),1);
        %for j = 1:length(indsIncluded)
        %    testPoint = [imX(indsIncluded(j)) imY(indsIncluded(j)) imZ(indsIncluded(j))];
        %    anglesToOurPoint(j) = acos(dot(intersectionPoint(ref,:,setC),testPoint)/(norm(intersectionPoint(ref,:,setC))*norm(testPoint))); 
        %end
        
        %from wikipedia...
        %1/sigma*sqrt(2*pi)*exo(-x^2/2/sigma^2)
        %fwhm = 2.355*sigma
        %acceptance angle is equivelent to FWHM
        sigma = acceptanceAngs(ref)/180*pi/2.355;
        gaussMod = 1/sigma/sqrt(2*pi)*exp(-anglesToOurPoint.^2/2/sigma^2); 

        if simBothEyes
           cordsToUse = voroniCordsBoth{i};
        else
           cordsToUse = voroniCordsLeft{i};
        end
        %does not account for sensitvity
        newImage(cordsToUse) = wmean(workingIm(indsIncluded),gaussMod);
        if sensitivityScale
            newImage(cordsToUse) = newImage(cordsToUse)*sensitivites(ref);
        end
        if sensitvityScaleByColSat
            satImg(cordsToUse) = sensitivites(ref);
        end
    end
    
    %build a hsv image and convert bak to rgb
    if sensitvityScaleByColSat
        newImgTemp = zeros(length(imAz), length(imEl), 3); 
        newImgTemp(:,:,1) = 0.16; %scale intensity in yellow...
        satImg = satImg-min(satImg(:));
        satImg = satImg./max(satImg(:)); %scale up to 0.5 saturation
        if invertSensColScale
            satImg = satImg*-1+1;
        end
        newImgTemp(:,:,2) = satImg;
        newImgTemp(:,:,3) = newImage/max(newImage(:));
        newImage = hsv2rgb(newImgTemp);
    end
    
    simulatedImF = figure; imshow_w_by_h(newImage/max(newImage(:))); hold on
     if plotBinoLine
        for i = 1:useBinoRegTo
            plot(binoImLines{i,1}, binoImLines{i,2}, 'c-', 'linewidth',2);
        end
    end
    
    if plotRightEyeLine
        for i = 1: useRightRegTo
            plot(leftImLines{i,1}, leftImLines{i,2}, 'g-', 'linewidth',2);
        end

        for i = 1:useRightRegTo
            plot(rightImLines{i,1}, rightImLines{i,2}, 'b-', 'linewidth',2);
        end
    end
    
    if displayTitles
       title(sprintf('Simulated insect view of image')); 
    end
    if displayAxisLabelsOnUnWrappedImages
       xlabel('Azimuth(^o)');
       ylabel('Elevation(^o)');
       axis on
       set(gca, 'Xtick', [1 imageHeight/2 imageHeight imageHeight*3/2 imageHeight*2], 'XTickLabel', {'\pm180', '-90', '0', '90', '\pm180'},...
           'Ytick', [1 imageHeight*0.25 imageHeight/2 imageHeight*0.75 imageHeight], 'YTickLabel', {'90', '45', '0', '-45', '-90'},'box','off');
    end
end

if saveData
   cd(saveDirectory);
   save(sprintf('%s.mat',saveFileName), 'lensDiamteres', 'interOAngle', 'eyeParam');
end

if saveImages
    cd(saveImDirectory);
%     saveas(vFOnSphereF, sprintf('%s_vFOnSphereF.pdf',saveFileName));
    
    saveas(plotLensDEyeF, sprintf('%s_plotLensDEyeF.tif',saveFileName));
    saveas(plotIOEyeF, sprintf('%s_plotIOEyeF.tif',saveFileName));
%     saveas(plotPEyeF, sprintf('%s_plotPEyeF.tif',saveFileName));
    
    saveas(diamUWFig, sprintf('%s_diamUWFig.svg',saveFileName));
    saveas(interOUWFig, sprintf('%s_interOUWFig.svg',saveFileName));
%     saveas(lensPUWFig, sprintf('%s_lensPUWFig.svg',saveFileName));
    
%     saveas(PColourbarFig, sprintf('%s_PColourbarFig.pdf',saveFileName));
    saveas(IOColourbarFig, sprintf('%s_IOColourbarFig.pdf',saveFileName));
    saveas(DColorbar, sprintf('%s_DColorbar.pdf',saveFileName));
    
    %may wish to change plotting commands for these
    saveas(linesOnOrigF, sprintf('%s_linesOnOrigF.svg',saveFileName));
    saveas(simulatedImF, sprintf('%s_simulatedImF.svg',saveFileName));
end