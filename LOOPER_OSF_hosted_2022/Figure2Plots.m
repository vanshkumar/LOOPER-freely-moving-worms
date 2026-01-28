
% run romoAnalysisPipeline_forFig3.m first
romoAnalysisPipeline_forFig3

if size(firingRates,1) > 1300
    firingRates = reshape(full(firingRates_sparse), firingRates_size);
    spikeRasters = reshape(full(spikeRasters_sparse), spikeRasters_size);
    
    % Neuron selection criteria used in the eLife paper
    D = size(trialNum,1);
    minN = min(reshape(trialNum(:,:,:), D, []), [], 2);
    meanFiringRate = mean(reshape(firingRatesAverage, D, []), 2);
    n = find(minN >= 5 & meanFiringRate < 50);
    
    firingRates = firingRates(n,:,:,:,:);
    firingRatesAverage = firingRatesAverage(n,:,:,:);
    spikeRasters = spikeRasters(n,:,:,:,:);
    trialNum = trialNum(n,:,:);
end


%%

% monkeySetupForFig2
load('monkeyDataForFigure3.mat');

F1Start = 0;
F2Start = 3.5;

trialStart = -2.5;
trialEnd = 5.5;

rasterTime = time;

plotNeurons = 1:size(firingRates,1);

plotRaster = squeeze(-spikeRasters(plotNeurons,1,1,:,5));
plotRaster = convn(plotRaster, [1 1 1], 'same');
plotRaster(plotRaster < -1) = -1;

figureHandle = figure(1);
figureHandle.Renderer='Painters';
clf;
hold on;
imagesc(time,1:length(plotNeurons),plotRaster);
colormap(gray);
xlim([trialStart, trialEnd])
ylim([1 length(plotNeurons)]);
yLimit = ylim;
plot([F1Start F1Start], [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F1Start F1Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F2Start F2Start], [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F2Start F2Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');

%%

figureHandle = figure(2);
figureHandle.Renderer='Painters';
clf;
hold on;
imagesc(time,1:length(plotNeurons),squeeze(firingRates(plotNeurons,1,1,:,5)));
colormap(parula);
xlim([trialStart, trialEnd])
ylim([1 length(plotNeurons)]);
yLimit = ylim;
plot([F1Start F1Start], [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F1Start F1Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F2Start F2Start], [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F2Start F2Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');

%%

% figureHandle = figure(2);
% figureHandle.Renderer='Painters';
% clf;
% hold on;
% imagesc(time,1:length(plotNeurons),squeeze(-firingRates(plotNeurons,1,1,:,1)));
% colormap(gray);
% xlim([trialStart, trialEnd])
% ylim([1 length(plotNeurons)]);
% yLimit = ylim;
% plot([F1Start F1Start], [yLimit(1) yLimit(2)], 'Color', 'k');
% plot([F1Start F1Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');
% plot([F2Start F2Start], [yLimit(1) yLimit(2)], 'Color', 'k');
% plot([F2Start F2Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');

%%

% run pfcRomo_dpcaPipleine.m first

X = (firingRatesAverage(:,:)');
times = 1:900;
doZScore = 0;
usePCA = 1;
numDelays = 0;
delayTime = 0;
figurePlot = 1;

[pcaBasis, ~, ~, ~, explained] = pca(X, 'NumComponents', 10);
Z = X * pcaBasis;
dataDim = size(firingRatesAverage);

Zfull = reshape(Z', [size(pcaBasis,2) dataDim(2:end)]);

colors = jet(7);
colors = colors([1 2 4 5 6 7],:);
colors(3,:) = colors(3,:)*0.6;
colors(4,:) = [238 210 2]/256;

stimMarg = find(strcmp(margNames, 'Stimulus'));
decisionMarg = find(strcmp(margNames, 'Decision'));
interactionMarg = find(strcmp(margNames, 'Interaction'));

topStim = find(whichMarg == stimMarg, 1);
topDecision = find(whichMarg == decisionMarg, 1);
topInteraction = find(whichMarg == interactionMarg, 1);

% plotComponents = [topStim, topDecision, topInteraction];
plotComponents = [topStim, topDecision];
% plotComponents = [topStim, 12, 18];
times = 1:900;
numDelays = 0;
delayTime = 0;
doZScore = 0;
usePCA = 0;

allTraces = [];
for f=1:size(Zfull,2)
    for d = 1:size(Zfull,3)
        thisTrace = squeeze(Zfull(plotComponents, f, d, times));        
        if doZScore
            thisTrace = zscore(thisTrace, [], 2);
        end
        thisTrace = delayEmbed(thisTrace, numDelays, delayTime);
        
        allTraces = [allTraces thisTrace];
    end
end

if size(allTraces,1) < 3
    pcaBasis = eye(size(allTraces,1),3);
else
    [pcaBasis, ~] = pca(allTraces', 'NumComponents', 6);
end

thisTrace = squeeze(Zfull(plotComponents, 1, 1, times));
if doZScore
    thisTrace = zscore(thisTrace, [], 2);
end
thisTrace = delayEmbed(thisTrace, numDelays, delayTime);
thisTrace = thisTrace'*pcaBasis;

figureHandle = figure(4);
figureHandle.Renderer='Painters';
clf;
hold on;

plot(time(times), thisTrace');
% colormap(parula);
xlim([trialStart, trialEnd])
% ylim([1 length(plotComponents)]);
yLimit = ylim;
plot([F1Start F1Start], [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F1Start F1Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F2Start F2Start], [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F2Start F2Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');


%% Apply dPCA

pfcRomo_dpcaPipelineForFig2


%% Get dPCA data

app.SavedData = saveData;

numTrial = max(app.SavedData.TrialData);
trialLength = size(app.SavedData.FinalStream,1) / numTrial;

matSize = size(allRawFiringRates);
matData = permute(allRawFiringRates, [1, 4, 2, 3, 5]);
allTrials = reshape(matData, [matSize(1), matSize(4), matSize(2) * matSize(3) * matSize(5)]);

inputData = permute(allInputs, [1, 4, 2, 3, 5]);
allInputs = reshape(inputData, [1, matSize(4), matSize(2) * matSize(3) * matSize(5)]);

X = rawStream;
times = 1:900;
doZScore = 0;
usePCA = 1;
numDelays = 0;
delayTime = 0;
figurePlot = 1;
decimateAmount = 10;

dPCAInput = firingRatesAverage(:,:)';

stimMarg = find(strcmp(margNames, 'Stimulus'));
decisionMarg = find(strcmp(margNames, 'Decision'));
interactionMarg = find(strcmp(margNames, 'Interaction'));

topStim = find(whichMarg == stimMarg, 2);
topDecision = find(whichMarg == decisionMarg, 2);
topInteraction = find(whichMarg == interactionMarg, 1);

plotComponents = [topStim  topDecision ];

rawTrials = [];
for i = 1:size(allTrials,1)
    for j = 1:size(allTrials,3)
        rawTrials(i,:,j) = decimate(squeeze(allTrials(i,:,j)),decimateAmount);
    end
end
delayedTrials = [];
for i = 1:size(allTrials,3)
    thisTrace = rawTrials(:,:,i)';
    
    Xcen = bsxfun(@minus, thisTrace, mean(dPCAInput));
    dPCAStream = Xcen * W;
    dPCAStream = dPCAStream(:,plotComponents)';
    
    delayedTrials(:,:,i) = delayEmbed(dPCAStream, 22, 1, 0, 1);
end
% totalDelay = app.SavedData.PreprocessData.DelayCount * app.SavedData.PreprocessData.DelayTime;
finalTrials = delayedTrials(:,end-1-trialLength:end-2,:);
rawStream = reshape(finalTrials, [size(finalTrials,1), size(finalTrials,2)*size(finalTrials,3)])';




% finalStream = Z(:,plotComponents);

%% Dispaly results

MIN_CLUSTER_SIZE = 4;

app.SavedData = saveData;

numTrial = max(app.SavedData.TrialData);
trialLength = size(app.SavedData.FinalStream,1) / numTrial;
times = 1:1000;

if max(times) > trialLength
    times = times(1):trialLength;
end

trialIndicies = repmat(times, [1, numTrial]);
trialIndicies = trialIndicies + kron((0:numTrial-1)*trialLength, ones(size(times)));

finalStream = zscore(rawStream, [], 1);
% finalStream = rawStream;
[trialPCABasis, ~] = pca(finalStream(trialIndicies,:), 'NumComponents',3);

% clf;
% h(1) = subplot(2,1,1)
% imagesc((finalStream*trialPCABasis)')
% 
% finalStream = app.SavedData.FinalStream;
% [trialPCABasis, ~] = pca(finalStream(trialIndicies,:), 'NumComponents',3);
% 
% h(2) = subplot(2,1,2)
% imagesc((finalStream*trialPCABasis)')
% linkaxes(h, 'x');
% xlim(66*12+[1 66])

loopStarts = (0:numTrial-1)*trialLength+times(1);

finalStream(loopStarts,:) = nan;

figureHandle = figure(7);
figureHandle.Renderer='Painters';
clf;
hold on;

colors = jet(7);
colors = colors([1 2 4 5 6 7],:);
colors(3,:) = colors(3,:)*0.6;
colors(4,:) = [238 210 2]/256;
% colors(4,:) = colors(3,:);

lineColors = [1 3 6 1 3 6];
lineColors = colors(lineColors,:);

% colors = jet(app.SavedData.BestLoopCount);
% h = plot3(finalStream(loopStarts+1,:)*trialPCABasis(:,1), finalStream(loopStarts+1,:)*trialPCABasis(:,2), finalStream(loopStarts+1,:)*trialPCABasis(:,3), 'o', 'LineWidth', 0.5);

% for i = [6, 7]%app.SavedData.BestLoopCount    
% for i = [4, 5]%app.SavedData.BestLoopCount    
% for i = [1,2,5,6]%app.SavedData.BestLoopCount    
% for i = [3,7,8,10]%app.SavedData.BestLoopCount    

for i = 1:numTrial
    thisTrialIDs = (i-1)*trialLength+times;
    conditionID = mode(mod(i-1,6)+1);
    
    h = plot3(finalStream(thisTrialIDs,:)*trialPCABasis(:,1), finalStream(thisTrialIDs,:)*trialPCABasis(:,2), finalStream(thisTrialIDs,:)*trialPCABasis(:,3), 'LineWidth', 0.5, 'Color', lineColors(conditionID,:));
    h.Color(4) = 0.2;
end

for i = 1:app.SavedData.BestLoopCount    
    thisLoopIDs = find(app.SavedData.BestLoopAssignments(:,1) == i);
    thisLoopClusters = app.SavedData.BestLoopAssignments(thisLoopIDs,2);
    
    allIDs = find(app.SavedData.BestStateMap(:,1) == i);
    trialIDs = floor((allIDs-1)/trialLength) + 1;
    conditionID = mode(mod(trialIDs-1,6)+1);
%     conditionID = mode(floor((trialIDs-1)/10)+1);
    
    badIDs = find(app.SavedData.BestStateMap(:,1) == i);
    badIDs = setdiff(1:size(finalStream,1), badIDs);
    
    plotStream = finalStream;
    plotStream(badIDs,:) = nan;
    
%     if PLOT_DELAY_EMBEDDED
%         h = plot3(plotStream(trialIndicies,:)*trialPCABasis(:,1), plotStream(trialIndicies,:)*trialPCABasis(:,2), plotStream(trialIndicies,:)*trialPCABasis(:,3), 'LineWidth', 0.5, 'Color', lineColors(conditionID,:));
%     else
%         h = plot(plotStream(trialIndicies,1), plotStream(trialIndicies,2), 'LineWidth', 0.5, 'Color', lineColors(conditionID,:));
%     end
%     h.Color(4) = 0.2;
    
    meanTimes = [];
    clusterLengths = [];
    clusterMeans = [];
    clusterSTDs = [];
    for j = 1:length(thisLoopClusters)
        thisIndices = find(app.SavedData.BestStateMap(:,1) == i & app.SavedData.BestStateMap(:,2) == thisLoopClusters(j));

        meanTimes(j) = mode(mod(thisIndices, trialLength)+1);
        clusterLengths(j) = length(thisIndices);
        
        clusterMeans(j,:) = nanmean(plotStream(thisIndices, :)*trialPCABasis(:,1:3), 1);
        clusterSTDs(j,:) = nanstd(plotStream(thisIndices, :)*trialPCABasis(:,1:3), [], 1);
    end
    
%     startPoints = find(app.SavedData.BestStateMap(loopStarts,1) == i);
%     bestStartPoint = mode(app.SavedData.BestStateMap(loopStarts(startPoints),2));
%     
%     
    thisLoopIDs = 1:length(thisLoopClusters);
    clusterOrder = 1:length(thisLoopClusters);
%     clusterOrder = app.SavedData.BestLoopAssignments(thisLoopIDs,2);
    
%     testTimes = meanTimes;
%     testTimes(clusterLengths < 5) = 1000;
    [~, startCluster] = min(meanTimes);
    
    sortedOrder = [startCluster:length(clusterOrder) 1:startCluster-1];
    thisLoopIDs = thisLoopIDs(sortedOrder);
    
%     for j = 1:length(clusterOrder)
%         thisIndices = find(app.SavedData.BestLoopAssignments(:,1) == i & app.SavedData.BestLoopAssignments(:,2) == clusterOrder(sortedOrder(j)));
%     end
    
    meanTimes = meanTimes(sortedOrder);
    
    badTimes = find(meanTimes <= min(times) | meanTimes >= max(times) | (clusterLengths < MIN_CLUSTER_SIZE));
    goodTimes = setdiff(sortedOrder, badTimes);
    
%     thisTrace = app.SavedData.BestEmission(thisLoopIDs(goodTimes), :);
    clusterSTDs = clusterSTDs(thisLoopIDs(goodTimes),:);
    clusterSTDs = filterData(clusterSTDs', 1, [], 1, 0)';
    thisTrace = clusterMeans(thisLoopIDs(goodTimes),:);
    thisTrace = filterData(thisTrace', 0.5, [], 1, 0)';
    
    tubeMean = thisTrace(:,1:3);
    projectedSTDs = clusterSTDs(:,1:3);
    
    projectedDerivative = diff(tubeMean);
    projectedDerivative = [projectedDerivative(1,:); projectedDerivative];
    
    orthogonalSTDs = zeros(size(projectedSTDs,1),1);
    for j = 1:size(projectedDerivative,1)
%         normDerivative = projectedDerivative(j,:)/norm(projectedDerivative(j,:));
%         derivativeNullspace = null(normDerivative)';
%         
%         changeOfBasis = [normDerivative; derivativeNullspace];
%         changedSTD = changeOfBasis * projectedSTDs(j,:)';
%         
%         orthogonalSTDs(j) = norm(changedSTD(2:end));
        
        orthogonalSTDs(j) = norm(projectedSTDs(j,:));
    end
    
    orthogonalSTDs = filterData(orthogonalSTDs, 0.5, [], 1, 0);
    
    if conditionID > 3
        lineStyle = ':';
    else
        lineStyle = '-';
    end
    
    h = plot3(tubeMean(:,1), tubeMean(:,2), tubeMean(:,3), lineStyle, 'LineWidth', 4, 'Color', lineColors(conditionID,:));
    
    [X,Y,Z,V] = tubeplot(tubeMean(:,1), tubeMean(:,2), tubeMean(:,3), orthogonalSTDs, ones(1,size(tubeMean,1)),10,[0 0 1]);

    t = surf(X,Y,Z,V, 'FaceColor', lineColors(conditionID,:));
    t.EdgeAlpha = 0.0;
    t.FaceAlpha = 0.1;
end

%% Reconstruction

decimateAmount = 10;

times = 1:900;
usF1s = [1 4 6];

app.SavedData = saveData;

numTrial = max(app.SavedData.TrialData);
trialLength = size(app.SavedData.FinalStream,1) / numTrial;
plotTimes = 1:500;

if max(plotTimes) > trialLength
    plotTimes = plotTimes(1):trialLength;
end



allTrials = firingRatesAverage(:,useF1s,:,times);

reconstructTrials = [];
for i = 1:size(allTrials,2)
    for j = 1:size(allTrials,3)
        trialData = [];
        
        for k = 1:size(allTrials,1)
            trialData(k,:) = decimate(squeeze(allTrials(k,i,j,:)),decimateAmount);
        end
        
        reconstructTrials(:,i,j,:) = trialData;
    end
end
reconstructTrials = reconstructTrials(:,:,:,end-1-trialLength:end-2);


matSize = size(allRawFiringRates);
matData = permute(allRawFiringRates, [1, 4, 2, 3, 5]);
allTrials = reshape(matData, [matSize(1), matSize(4), matSize(2) * matSize(3) * matSize(5)]);

inputData = permute(allInputs, [1, 4, 2, 3, 5]);
allInputs = reshape(inputData, [1, matSize(4), matSize(2) * matSize(3) * matSize(5)]);

rawTrials = [];
for i = 1:size(allTrials,1)
    for j = 1:size(allTrials,3)
        rawTrials(i,:,j) = decimate(squeeze(allTrials(i,:,j)),decimateAmount);
    end
end
rawTrials = rawTrials(:,end-1-trialLength:end-2,:);
rawStream = reshape(rawTrials, [size(rawTrials,1), size(rawTrials,2)*size(rawTrials,3)])';


trialIndicies = repmat(plotTimes, [1, numTrial]);
trialIndicies = trialIndicies + kron((0:numTrial-1)*trialLength, ones(size(plotTimes)));

loopConditions = [];
for i = 1:app.SavedData.BestLoopCount    
    allIDs = find(app.SavedData.BestStateMap(:,1) == i);
    trialIDs = floor((allIDs-1)/trialLength) + 1;
    
    loopConditions(i) = mode(mod(trialIDs-1,6)+1);
end

clusterMeans = [];
clusterSTD = [];
clusterLengths = [];
meanTimes = [];
stateLookup = [];
for i = 1:size(app.SavedData.BestLoopAssignments,1)
    loopID = app.SavedData.BestLoopAssignments(i,1);
    phaseID = app.SavedData.BestLoopAssignments(i,2);
    
    thisIndices = find(app.SavedData.BestStateMap(:,1) == loopID & app.SavedData.BestStateMap(:,2) == phaseID);
    
    trialTimes = mod(thisIndices, trialLength)+1;
    meanTimes(i,:) = mode(trialTimes);
    clusterLengths(i,:) = length(thisIndices);

    clusterMeans(i,:) = nanmean(rawStream(thisIndices, :), 1);
    clusterSTD(i,:) = nanstd(rawStream(thisIndices, :), [], 1);
    
    stateLookup(loopID, phaseID) = i;
end

reconstructedStream = [];
conditionCounts = zeros(size(reconstructTrials,2), size(reconstructTrials,3));
for i = 1:size(rawTrials,3)
    conditionID = mod((i-1),6)+1;
    
    F1ID = mod(conditionID - 1,3)+1;
    F2ID = floor((conditionID-1)/3)+1;
    
    thisIDs = (i-1)*trialLength + (1:trialLength);
    thisLoopIDs = app.SavedData.BestStateMap(thisIDs,1);
    thisPhaseIDs = app.SavedData.BestStateMap(thisIDs,2);
    
    thisStateIDs = [];
    for j = 1:length(thisLoopIDs)
        thisStateIDs(j) = stateLookup(thisLoopIDs(j), thisPhaseIDs(j));
    end
    
    conditionCounts(F1ID, F2ID) = conditionCounts(F1ID, F2ID) + 1;
    
    reconstructedStream(:, F1ID, F2ID, :, conditionCounts(F1ID, F2ID)) = clusterMeans(thisStateIDs,:)';
end

Rsquareds = [];
for i = 1:size(reconstructedStream,2)
    for j = 1:size(reconstructedStream,3)
        for k = 1:size(reconstructedStream,1)
            for trialNumber = 1:size(reconstructedStream,5)
                thisReconstruction = squeeze(reconstructedStream(k,i,j,:,trialNumber));
                targetData = squeeze(reconstructTrials(k,i,j,:));
                
                Rsquareds(k,i,j,trialNumber) = corr(thisReconstruction, targetData);
            end
        end
    end
end

meanRsquared = nanmean(Rsquareds(:))

trialRsquareds = [];
for i = 1:size(reconstructedStream,2)
    for j = 1:size(reconstructedStream,3)
        for trialNumber = 1:size(reconstructedStream,5)
            thisReconstruction = squeeze(reconstructedStream(:,i,j,:,trialNumber));
            targetData = squeeze(reconstructTrials(:,i,j,:));

            trialRsquareds(i,j,trialNumber) = corr(thisReconstruction(:), targetData(:));
        end
    end
end

meanTrialRsquared = nanmean(trialRsquareds(:))

%% Plot traces

plotF1 = 3;
plotF2 = 2;
plotTrial = 5;

originalTime = time;
originalTime = decimate(originalTime(times),decimateAmount);
originalTime = originalTime(end-1-trialLength:end-2);

thisData = squeeze(reconstructTrials(plotNeurons,plotF1,plotF2,:));

distances = pdist(thisData, 'correlation');
distances(isnan(distances)) = 0;
heirarchy = linkage(distances, 'average');

leafOrder = optimalleaforder(heirarchy, distances);


figureHandle = figure(8);
figureHandle.Renderer='Painters';
clf;
hold on;
imagesc(originalTime,1:length(plotNeurons),squeeze(reconstructedStream(leafOrder(plotNeurons),plotF1,plotF2,:,plotTrial)));
colormap(parula);
xlim([min(originalTime), trialEnd])
ylim([1 length(plotNeurons)]);
yLimit = ylim;
plot([F1Start F1Start], [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F1Start F1Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F2Start F2Start], [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F2Start F2Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');
title('Reconstructed Data');

figureHandle = figure(9);
figureHandle.Renderer='Painters';
clf;
hold on;
imagesc(originalTime,1:length(plotNeurons),squeeze(reconstructTrials(leafOrder(plotNeurons),plotF1,plotF2,:)));
colormap(parula);
xlim([min(originalTime), trialEnd])
ylim([1 length(plotNeurons)]);
yLimit = ylim;
plot([F1Start F1Start], [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F1Start F1Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F2Start F2Start], [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F2Start F2Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');
title('Average Data');


plotRaster = squeeze(-spikeRasters(leafOrder(plotNeurons),plotF1,plotF2,:,plotTrial));
plotRaster = convn(plotRaster, [1 1 1], 'same');
plotRaster(plotRaster < -1) = -1;

figureHandle = figure(10);
figureHandle.Renderer='Painters';
clf;
hold on;
imagesc(time,1:length(plotNeurons),plotRaster);
colormap(gray);
xlim([min(originalTime), trialEnd])
ylim([1 length(plotNeurons)]);
yLimit = ylim;
plot([F1Start F1Start], [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F1Start F1Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F2Start F2Start], [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F2Start F2Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');
