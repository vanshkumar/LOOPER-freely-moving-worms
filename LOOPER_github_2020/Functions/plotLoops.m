%% Display results

PLOT_DELAY_EMBEDDED = 1;

PLOT_PHASE = 0;

numTrial = max(saveData.TrialData);
reducedSize = size(saveData.TrialData,2) - size(saveData.FinalStream,1);
reducedSize = reducedSize / numTrial;
% trialLength = size(saveData.FinalStream,1) / numTrial;
% times = 1:trialLength;

% if max(times) > trialLength
%     times = times(1):trialLength;
% end

% trialIndicies = repmat(times, [1, numTrial]);
% trialIndicies = trialIndicies + kron((0:numTrial-1)*trialLength, ones(size(times)));

finalStream = saveData.FinalStream;

if size(saveData.ClusterMeans,2) > 2
    [trialPCABasis, ~] = pca(saveData.ClusterMeans, 'NumComponents',3);
else
    trialPCABasis = eye(size(saveData.ClusterMeans,2));
end
% trialPCABasis = eye(size(finalStream,2), 3);

loopStarts = [0, saveData.RawTrialSwitches] - [0:numTrial-1]*reducedSize + 1;

longestTrial = max(diff([loopStarts size(finalStream,2)]));

% times = 1:;
% 
% if max(times) > trialLength
%     times = times(1):trialLength;
% end
% 
% trialIndicies = loopStarts;
% trialIndicies = trialIndicies + kron((0:numTrial-1)*trialLength, ones(size(times)));

finalStream(loopStarts,:) = nan;

figureHandle = figure(10000);
figureHandle.Renderer='Painters';
clf;
hold on;
title('LOOPER loops (tubes in model space)');

tempLines = lines;
clear lines
if saveData.BestLoopCount < 8
    colors = lines(saveData.BestLoopCount);
else
    colors = jet(saveData.BestLoopCount);
end
lines = tempLines;
%colors = colors([1 2 4 5 6 7],:);
%colors(3,:) = colors(3,:)*0.6;
%colors(4,:) = [238 210 2]/256;
% colors(4,:) = colors(3,:);

lineColors = 1:saveData.BestLoopCount;
lineColors = colors(lineColors,:);

allClusterMeans = [];
for i = 1:saveData.BestLoopCount
    thisLoopIDs = find(saveData.BestLoopAssignments(:,1) == i);
    thisLoopClusters = saveData.BestLoopAssignments(thisLoopIDs,2);
    
    allIDs = find(saveData.BestStateMap(:,1) == i);
%     trialIDs = floor((allIDs-1)/trialLength) + 1;
    conditionID = i;
    
    badIDs = find(saveData.BestStateMap(:,1) == i);
    badIDs = setdiff(1:size(finalStream,1), badIDs);
    
    plotStream = finalStream;
    plotStream(badIDs,:) = nan;
    
    if PLOT_DELAY_EMBEDDED
        if size(trialPCABasis,2) > 2
            h = plot3(plotStream(:,:)*trialPCABasis(:,1), plotStream(:,:)*trialPCABasis(:,2), plotStream(:,:)*trialPCABasis(:,3), 'LineWidth', 1, 'Color', lineColors(conditionID,:));
        else
            h = plot(plotStream(:,:)*trialPCABasis(:,1), plotStream(:,:)*trialPCABasis(:,2), 'LineWidth', 1, 'Color', lineColors(conditionID,:));
        end
    else
        h = plot3(plotStream(:,1), plotStream(:,2), plotStream(:,3), 'LineWidth', 1, 'Color', lineColors(conditionID,:));
    end
    h.Color(4) = 0.3;
    
    meanTimes = [];
    clusterLengths = [];
    clusterMeans = [];
    clusterSTDs = [];
    for j = 1:length(thisLoopClusters)
        thisIndices = find(saveData.BestStateMap(:,1) == i & saveData.BestStateMap(:,2) == thisLoopClusters(j));

        loopEnds = [loopStarts, size(saveData.FinalStream, 1)];
        
        thisTimes = thisIndices;
        currentIndex = 1;
        oldIndex = 1;
        for k = 2:length(loopEnds)
            validTimes = thisTimes + 1 - loopEnds(k);
            currentIndex = currentIndex + find(validTimes(currentIndex:end) > 0, 1) - 2;
            if k == length(loopEnds)
                currentIndex = length(thisTimes);
            end
            if currentIndex == 0
                currentIndex = 1;
                continue;
            end
            
            thisTimes(oldIndex:currentIndex) = thisTimes(oldIndex:currentIndex) - loopEnds(k-1);
            
            oldIndex = currentIndex+1;
        end
        
        meanTimes(j) = mode(thisTimes);
        clusterLengths(j) = length(thisIndices);
        
        clusterMeans(j,:) = nanmean(plotStream(thisIndices, :), 1);
        clusterSTDs(j,:) = nanstd(plotStream(thisIndices, :), [], 1);
        
        allClusterMeans(i, thisLoopClusters(j),:) = clusterMeans(j,:);
    end
    
    thisLoopIDs = 1:length(thisLoopClusters);
    clusterOrder = 1:length(thisLoopClusters);
    
    [~, startCluster] = min(meanTimes);
    
    sortedOrder = [startCluster:length(clusterOrder) 1:startCluster-1];
    thisLoopIDs = thisLoopIDs(sortedOrder);
    
    meanTimes = meanTimes(sortedOrder);
    
    badTimes = [];%find(meanTimes <= min(times) | meanTimes >= max(times));
    goodTimes = setdiff(sortedOrder, badTimes);
    
    clusterSTDs = clusterSTDs(thisLoopIDs(goodTimes),:);
    thisTrace = clusterMeans(thisLoopIDs(goodTimes),:);
    
    if length(goodTimes) <= 1
        continue;
    end
    
    thisTrace = filterData(thisTrace', 0.5, [], 1, 0)';
    
    if PLOT_DELAY_EMBEDDED
        maxDim = min(3, size(trialPCABasis,2));
        
        tubeMean = thisTrace*trialPCABasis(:,1:maxDim);
        projectedSTDs = clusterSTDs*trialPCABasis(:,1:maxDim);
    else
        tubeMean = [thisTrace(:,1:3) zeros(size(thisTrace,1), 1)];
        projectedSTDs = [clusterSTDs(:,1:3) zeros(size(thisTrace,1), 1)];
    end
    
    projectedDerivative = diff(tubeMean);
    projectedDerivative = [projectedDerivative(1,:); projectedDerivative];
    
    orthogonalSTDs = zeros(size(projectedSTDs,1),1);
    for j = 1:size(projectedDerivative,1)
        orthogonalSTDs(j) = norm(projectedSTDs(j,:));
    end
    
    orthogonalSTDs = filterData(orthogonalSTDs, 0.5, [], 1, 2);
    
    if size(trialPCABasis,2) > 2
        [X,Y,Z,V] = tubeplot(tubeMean(:,1), tubeMean(:,2), tubeMean(:,3), orthogonalSTDs, ones(1,size(tubeMean,1)),10,[0 0 1]);
    else
        [X,Y,Z,V] = tubeplot(tubeMean(:,1), tubeMean(:,2), tubeMean(:,1)*0, orthogonalSTDs, ones(1,size(tubeMean,1)),10,[0 0 1]);
    end
    
    if PLOT_PHASE    
        colors = hsv(256);
        discreteTimes = floor(meanTimes / longestTrial * 256) + 1;
        colors = reshape(colors(discreteTimes,:), [size(X,1) 1 3]);
        colors = repmat(colors, [1 size(X,2), 1]);
    
        t = surf(X,Y,Z,colors);
        t.EdgeAlpha = 0.2;
        t.FaceAlpha = 0.5;
    else        
        t = surf(X,Y,Z,V, 'FaceColor', lineColors(conditionID,:));
        t.EdgeColor = lineColors(conditionID,:);
        t.EdgeAlpha = 0.5;
        t.FaceAlpha = 0.2;
    end
    
    
end
