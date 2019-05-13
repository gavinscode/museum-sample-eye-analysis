function [ voroniLabels ] = assignSubsToVoroniDiagram(xyz, sX, sY, sZ )
    %Written by Gavin Taylor, 2017. MIT License

    voroniLabels = zeros(length(sX),1); 

    xyz = bsxfun(@rdivide, xyz, sqrt(sum(xyz.^2,1)));
    [~, ~, voronoiboundary] = voronoisphere(xyz);
    
    for i = 1:size(voronoiboundary,1)
        X = voronoiboundary{i};
        closerPts = X*0.9;
        furtherPts = X*1.1;
        fullPts = [closerPts furtherPts]';

        if sum(isnan(fullPts(:))) | sum(isinf(fullPts(:)))
            indsToR = find( isnan(fullPts(:,1)) | isinf(fullPts(:,1)));
            fullPts(indsToR,:) = [];
        end
        bFace = boundary(fullPts(:,1), fullPts(:,2), fullPts(:,3),0.25);
        xR = [min(fullPts(:,1)) max(fullPts(:,1))];
        yR = [min(fullPts(:,2)) max(fullPts(:,2))];
        zR = [min(fullPts(:,3)) max(fullPts(:,3))];
        toTest = find(sX >= xR(1) & sX <= xR(2) & sY >= yR(1) & sY <= yR(2) & sZ >= zR(1) & sZ <= zR(2));
        cellFOV = find(intriangulation(fullPts, bFace, [sX(toTest) sY(toTest) sZ(toTest)]));
        
        voroniLabels(toTest(cellFOV)) = i; 

%         figure; hold on
%         plot3(fullPts(:,1), fullPts(:,2), fullPts(:,3),'bx')
%         plot3(sX(toTest), sY(toTest), sZ(toTest),'go')
%         plot3(sX(:), sY(:), sZ(:),'r.')
%         patch(struct('faces', bFace,'vertices',fullPts));
    end
    
%     f=figure;
%     set(f,'Renderer','zbuffer');
%     ax = axes('Parent', f);
%     hold(ax, 'on');
%     axis(ax,'equal');
% 
%     plot3(ax, xyz(1,:),xyz(2,:),xyz(3,:),'wo');
%     clmap = cool();
%     ncl = size(clmap,1);
%     for k = 1:size(voronoiboundary,1)
%         X = voronoiboundary{k};
%         cl = clmap(mod(k,ncl)+1,:);
%         fill3(X(1,:),X(2,:),X(3,:),cl,'Parent',ax,'EdgeColor','w');
%     end
%     axis(ax,'equal');
end
%to plot full veroni lines on image
% for i = 1:size(voronoiboundary,1)
%     X = voronoiboundary{i};
%     tempEl = zeros(size(X,2),1);
%     tempAz = zeros(size(X,2),1);
% 
%     vorEls = acos(X(2,:)')/pi*180-90;
%     vorAzs = atan2(X(1,:)',-X(3,:)')/pi*180;
% 
%     for j = 1:size(X,2)
%         [~, tempAz(j)] = min(abs(imAz-vorAzs(j)));
%         [~, tempEl(j)] = min(abs(imEl-vorEls(j)));
%     end
%     plot(tempAz, tempEl, 'm-')
% end
