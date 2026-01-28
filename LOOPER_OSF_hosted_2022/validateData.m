MINIMUM_STATE_TIME = 0 *(2+1);

app.SavedData = saveData;

numTrial = max(saveData.TrialData);
reducedSize = size(saveData.TrialData,2) - size(saveData.FinalStream,1);
reducedSize = reducedSize / numTrial;

finalStream = saveData.FinalStream;

if size(saveData.ClusterMeans,2) > 2
    [trialPCABasis, ~] = pca(saveData.ClusterMeans, 'NumComponents',3);
else
    trialPCABasis = eye(size(saveData.ClusterMeans,2));
end

loopStarts = [0, saveData.RawTrialSwitches] - [0:numTrial-1]*reducedSize + 1;

longestTrial = max(diff([loopStarts size(finalStream,2)]));

finalStream(loopStarts,:) = nan;

loopStarts(end+1) = size(finalStream,1)+1;

if SHOULD_VALIDATE    
    rawData = convertToCell(allTrials);

    lastEnd = 0;
    rawTrialData = [];
    for i = 1:length(rawData)
        rawTrialData(lastEnd + (1:size(rawData{i}, 2))) = i;

        lastEnd = lastEnd + size(rawData{i}, 2);
    end

    rawData = mergeData(rawData);

    [tempData, trialData, procesedTrialSwitches] = preprocessData(rawData, [], 0, [], 0, rawTrialData, 0, true, saveData.PreprocessData.Smoothing, saveData.PreprocessData.ZScore, saveData.PreprocessData.DelayTime, saveData.PreprocessData.DelayCount, saveData.DataMean, saveData.DataSTD);
    badIndicies = [procesedTrialSwitches procesedTrialSwitches+1];
    badIndicies(badIndicies > size(tempData,2)) = [];
    tempData(:,badIndicies) = [];


    clusterDistances = [];
    finalStream = app.SavedData.FinalStream;
    for i = 1:size(app.SavedData.BestEmission,1)
        thisLoopPosition = app.SavedData.BestLoopAssignments(i,:);
        thisIndices = find(ismember(app.SavedData.BestStateMap, thisLoopPosition, 'rows'));

        stds = std(finalStream(thisIndices,:), [], 1);

        thisEmission = mean(finalStream(thisIndices,:), 1);

        clusterStream = (thisEmission - tempData') ./ stds;           

        clusterDistances(i,:) = sum(clusterStream.^2,2);
    end

    [~, bestClusters] = min(clusterDistances, [], 1);

    loopIDs = app.SavedData.BestLoopAssignments(bestClusters,1);
    phaseIDs = app.SavedData.BestLoopAssignments(bestClusters,2);

    numTrial = length(allTrials);

    MINIMUM_STATE_TIME = 5;
else    
    loopIDs = app.SavedData.BestStateMap(:,1);
    phaseIDs = app.SavedData.BestStateMap(:,2);
end
if MINIMUM_STATE_TIME > 0
    loopIDs = colfilt(loopIDs, [MINIMUM_STATE_TIME 1], 'sliding', @mode);
    phaseIDs = colfilt(phaseIDs, [MINIMUM_STATE_TIME 1], 'sliding', @mode);
end

if ~exist('loopOrder', 'var') || isempty(loopOrder)
    loopOrder = 1:max(loopIDs);
end

clear lines;
lineColors = lines(length(loopOrder));

% trialLoopIDs = [];
% trialConditionIDs = [];

allPhases = [];

loopAssignments = nan(numTrial, loopStarts(2)-loopStarts(1));
phaseAssignments = nan(numTrial, loopStarts(2)-loopStarts(1));

figureHandle = figure(1);
figureHandle.Renderer='Painters';
clf;
hold on;
title('Loop ID over time (by trial)');
legendLines = [];
for i = 1:numTrial    
    conditionID = floor((i-1)/totalTrials)+1;
    trialIDs = loopStarts(i):loopStarts(i+1)-1;

    IDs = loopIDs(trialIDs);
    IDs(IDs == 0) = length(loopOrder);

    plotTime = 1:length(IDs);

    h = plot(plotTime, loopOrder(IDs) + normrnd(0,0.03,size(trialIDs)), 'Color', lineColors(conditionID,:));
    h.Color(4) = 0.2;
%     scatter(plotTime, loopOrder(IDs) + normrnd(0,0.03,size(trialIDs)), 32, phaseIDs(trialIDs));
    legendLines(conditionID) = h;

    allPhases(i,:) = phaseIDs(trialIDs);
    
    loopAssignments(i, :) = IDs;

%     trialLoopIDs(:,i) = loopOrder(IDs);
%     trialConditionIDs(i) = conditionID;

%     transitionCheck = colfilt(loopOrder(IDs), [1 2], 'sliding', @mean);
%     transitionCheck(transitionCheck == loopOrder(IDs)) = 0;
%     transitionPoints = find(transitionCheck);
%     transitionPoints = [transitionPoints, transitionPoints+1];
%     transitionPoints(transitionPoints > trialLength) = [];
%     
%     transitionData = nan(size(loopOrder(IDs)));
%     transitionData(transitionPoints) = loopOrder(IDs(transitionPoints));

%     h = plot(plotTime, transitionData + normrnd(0,0.03,size(trialIDs)), lineColors{conditionID});
%     h.Color(4) = 0.2;
end
clear colormap;
colormap(hsv(256));

phaseAssignments = allPhases;
