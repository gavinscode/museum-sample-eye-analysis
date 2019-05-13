%Written by Gavin Taylor, 2017. MIT License

%Run to create group histogram plots, and also to set colour limit bounds
%for main script.

clear; close all; clc

cd('Insert data directory here');
dataToGet = {'Insert name of data files here', 'One cell per file'};

%%%set desired number of color bins here and set this parameter to []
%in main to use same color map for all
numColBins = 10; %actually +1

%for plotting 
sFac = 4; disPara = 0.5;
loadedData = cell(length(dataToGet),1);

for i = 1:length(dataToGet)
   loadedData{i} = load(dataToGet{i}); 
end

% get means,stds and ranges for each
for i = 1:length(dataToGet)
    tempDat = loadedData{i};
    fprintf('\n%s\n',dataToGet{i});
    fprintf('Facet D mean: %.1f std: %.1f min: %.1f max: %.1f\n', mean(tempDat.lensDiamteres), std(tempDat.lensDiamteres), min(tempDat.lensDiamteres), max(tempDat.lensDiamteres));
    fprintf('Inter O mean: %.1f std: %.1f min: %.1f max: %.1f\n', mean(tempDat.interOAngle), std(tempDat.interOAngle), min(tempDat.interOAngle), max(tempDat.interOAngle));
    fprintf('Eye P mean: %.1f std: %.1f min: %.1f max: %.1f\n', mean(tempDat.eyeParam), std(tempDat.eyeParam), min(tempDat.eyeParam), max(tempDat.eyeParam));
end
   
% plot histograms
close all
interOHistFig = figure; hold on;
diameterHistFig = figure; hold on
eyeParamHistFig = figure; hold on
legendDets = cell(length(dataToGet),1);
for i = 1:length(dataToGet)
    tempDat = loadedData{i};
    
    figure(interOHistFig);
    [n, x] = hist(tempDat.interOAngle,0:disPara:50);
    plot(x,smooth(n,sFac))

    figure(diameterHistFig);
    [n, x] = hist(tempDat.lensDiamteres,0:disPara:30);
    plot(x,smooth(n,sFac))

    figure(eyeParamHistFig);
    [n, x] = hist(tempDat.eyeParam,0:disPara:30);
    plot(x,smooth(n,sFac))

    tmpName = dataToGet{i};
    legendDets{i} = sprintf('%s', tmpName(1:end-4));
end

figure(interOHistFig);
ylim([0 80]); xlim([0 30]);
xlabel('Interommatidial angle (^o)'); ylabel('# facets');
set(gca, 'TickDir','out', 'XTick', [0 10 20 30], 'YTick', 0:20:80);
pos = get(gcf, 'Position');
pos(1:2) = 100; pos(3) = 500*1.62; pos(4) = 500;
% set(gcf, 'Position', pos)
set(gcf,'Position',[50 50 1175 365])

figure(diameterHistFig);
ylim([0 80]); xlim([0 30]);
xlabel('Facet size (um)'); ylabel('# facets');
set(gca, 'TickDir','out', 'XTick', [0 10 20 30], 'YTick', 0:20:80);
pos = get(gcf, 'Position');
pos(1:2) = 100; pos(3) = 500*1.62; pos(4) = 500;
% set(gcf, 'Position', pos)
set(gcf,'Position',[50 50 1175 365])
legend(legendDets);
     
figure(eyeParamHistFig);
ylim([0 80]); xlim([0 15]);
xlabel('Eye parameter (um.rad)'); ylabel('# facets');
set(gca, 'TickDir','out', 'XTick', [0 5 10 15], 'YTick', 0:20:80);
pos = get(gcf, 'Position');
pos(1:2) = 100; pos(3) = 500*1.62; pos(4) = 500;
% set(gcf, 'Position', pos)
set(gcf,'Position',[50 50 1175 365])

% set colour map bounds over range of eyes
fullDList = [];
fullIOList = [];
fullPList = [];

for i = 1:length(dataToGet)
    tempDat = loadedData{i};
    fullDList = [fullDList tempDat.lensDiamteres'];
    fullIOList = [fullIOList tempDat.interOAngle'];
    fullPList = [fullPList tempDat.eyeParam'];
end

 [n, x] = hist(fullDList,100);
 cumSAngle = cumsum(n)/sum(n);
 %do proportional division
 xDLower = zeros(numColBins+1,1);
 xDLower(1) = x(1)-1; xDLower(end) = x(end)+1;
 for i = 1:numColBins-1
     temp = find(cumSAngle > i/numColBins);
     xDLower(i+1) = x(temp(1));
 end
 %do linear divsion
temp = find(cumSAngle > 0.01); lower = x(temp(1));
temp = find(cumSAngle > 0.99); upper = x(temp(1));
xDLower = lower:(upper-lower)/(numColBins):upper

 [n, x] = hist(fullIOList,100);
 cumSAngle = cumsum(n)/sum(n);
  %do proportional division
 xIOLower = zeros(numColBins+1,1);
 xIOLower(1) = x(1)-1; xIOLower(end) = x(end)+1;
 for i = 1:numColBins-1
     temp = find(cumSAngle > i/numColBins);
     xIOLower(i+1) = x(temp(1));
 end
 %do linear divsion
temp = find(cumSAngle > 0.01); lower = x(temp(1));
temp = find(cumSAngle > 0.99); upper = x(temp(1));
xIOLower = lower:(upper-lower)/(numColBins):upper

 [n, x] = hist(fullPList,100);
 cumSAngle = cumsum(n)/sum(n);
  %do proportional division
 xPLower = zeros(numColBins+1,1);
 xPLower(1) = x(1)-1; xPLower(end) = x(end)+1;
 for i = 1:numColBins-1
     temp = find(cumSAngle > i/numColBins);
     xPLower(i+1) = x(temp(1));
 end
 %do linear divsion
temp = find(cumSAngle > 0.01); lower = x(temp(1));
temp = find(cumSAngle > 0.99); upper = x(temp(1));
xPLower = lower:(upper-lower)/(numColBins):upper
