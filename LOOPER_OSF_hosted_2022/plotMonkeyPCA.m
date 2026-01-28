%% Display results

load('monkey3F1.mat');

PLOT_DELAY_EMBEDDED = 1;
MIN_CLUSTER_SIZE = 4;

app.SavedData = saveData;

% useData = app.SavedData.FinalStream;
useData = app.SavedData.RawData';

numTrial = max(app.SavedData.TrialData);
trialLength = size(useData,1) / numTrial;
times = 1:1000;

if max(times) > trialLength
    times = times(1):trialLength;
end

trialIndicies = repmat(times, [1, numTrial]);
trialIndicies = trialIndicies + kron((0:numTrial-1)*trialLength, ones(size(times)));

finalStream = useData;

loopStarts = (0:numTrial-1)*trialLength+times(1);

finalStream(loopStarts,:) = nan;

meanData = zeros(trialLength, size(finalStream,2));
for i = 1:numTrial
    thisTrialIDs = (i-1)*trialLength+times;
    
    meanData = meanData + finalStream(thisTrialIDs,:);
end
meanData = meanData / numTrial;

for i = 1:numTrial
    thisTrialIDs = (i-1)*trialLength+times;
    
    finalStream(thisTrialIDs,:) = finalStream(thisTrialIDs,:) - meanData;
end

[trialPCABasis, ~, ~, ~,explained] = pca(finalStream(trialIndicies,:), 'NumComponents',3);
% trialPCABAsis = pcaBasis;

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
    
    style = '-';
    if conditionID > 3
        style = ':';
    end
    
    if PLOT_DELAY_EMBEDDED
        h = plot3(finalStream(thisTrialIDs,:)*trialPCABasis(:,1), finalStream(thisTrialIDs,:)*trialPCABasis(:,2), finalStream(thisTrialIDs,:)*trialPCABasis(:,3), style, 'LineWidth', 1, 'Color', lineColors(conditionID,:));
    else
        h = plot(finalStream(thisTrialIDs,1), finalStream(thisTrialIDs,2), style, 'LineWidth', 1, 'Color', lineColors(conditionID,:));
    end
    h.Color(4) = 0.6;
end

%%

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
        
        clusterMeans(j,:) = nanmean(plotStream(thisIndices, :), 1);
        clusterSTDs(j,:) = nanstd(plotStream(thisIndices, :), [], 1);
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
    
    if PLOT_DELAY_EMBEDDED
%         plot3(thisTrace*trialPCABasis(:,1), thisTrace*trialPCABasis(:,2), thisTrace*trialPCABasis(:,3), 'LineWidth', 2, 'Color', colors(i,:))
    
        tubeMean = thisTrace*trialPCABasis(:,1:3);
        projectedSTDs = clusterSTDs*trialPCABasis(:,1:3);
    else
%         plot(thisTrace(:,1), thisTrace(:,2), 'LineWidth', 2, 'Color', colors(i,:))
        
        tubeMean = [thisTrace(:,1:2) zeros(size(thisTrace,1), 1)];
        projectedSTDs = [clusterSTDs(:,1:2) zeros(size(thisTrace,1), 1)];
    end
    
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
%     h.Color(4) = 0.5;
    
    [X,Y,Z,V] = tubeplot(tubeMean(:,1), tubeMean(:,2), tubeMean(:,3), orthogonalSTDs, ones(1,size(tubeMean,1)),10,[0 0 1]);

    t = surf(X,Y,Z,V, 'FaceColor', lineColors(conditionID,:));
    t.EdgeAlpha = 0.0;
    t.FaceAlpha = 0.1;
end

