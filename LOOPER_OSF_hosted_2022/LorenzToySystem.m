
%% Simulate lorenz

noise = 1;
dt = 0.01;

sigma = 10;
beta = 8/3;
rho = 28;

startingPoint = [2, 3, 4];
traces = startingPoint;

for t = 2:10000
    dx = sigma*(traces(t-1,2) - traces(t-1,1));
    dy = traces(t-1,1)*(rho - traces(t-1,3)) - traces(t-1,2);
    dz = traces(t-1,1)*traces(t-1,2) - beta*traces(t-1,3);
    
    traces(t,1) = traces(t-1,1) + dx * dt;
    traces(t,2) = traces(t-1,2) + dy * dt;
    traces(t,3) = traces(t-1,3) + dz * dt;
    
    traces(t,:) = traces(t,:) + normrnd(0, noise, size(startingPoint));
end


finalTrials = [];
for i = 1:size(traces,2)
    finalTrials(i,:) = decimate(traces(:,i), 4);
end

figure(1);
clf;
h = plot3(finalTrials(1,:), finalTrials(2,:), finalTrials(3,:));
h.Color(4) = 0.1;

%% Display results

PLOT_DELAY_EMBEDDED = 1;

app.SavedData = saveData;

numTrial = max(app.SavedData.TrialData);
trialLength = size(app.SavedData.FinalStream,1) / numTrial;
times = 1:trialLength;

if max(times) > trialLength
    times = times(1):trialLength;
end

trialIndicies = repmat(times, [1, numTrial]);
trialIndicies = trialIndicies + kron((0:numTrial-1)*trialLength, ones(size(times)));

finalStream = app.SavedData.FinalStream;

% [trialPCABasis, ~] = pca(finalStream(trialIndicies,:), 'NumComponents',3);
trialPCABasis = eye(size(finalStream,2), 3);

loopStarts = (0:numTrial-1)*trialLength+1;

finalStream(loopStarts,:) = nan;

figureHandle = figure(1);
figureHandle.Renderer='Painters';
clf;
hold on;

colors = jet(7);
colors = colors([1 2 4 5 6 7],:);
colors(3,:) = colors(3,:)*0.6;
colors(4,:) = [238 210 2]/256;
% colors(4,:) = colors(3,:);

lineColors = 1:app.SavedData.BestLoopCount;
lineColors = colors(lineColors,:);

% colors = jet(app.SavedData.BestLoopCount);
% h = plot3(finalStream(loopStarts+1,:)*trialPCABasis(:,1), finalStream(loopStarts+1,:)*trialPCABasis(:,2), finalStream(loopStarts+1,:)*trialPCABasis(:,3), 'o', 'LineWidth', 0.5);

% for i = [6, 7]%app.SavedData.BestLoopCount    
% for i = [4, 5]%app.SavedData.BestLoopCount    
% for i = [1,2,5,6]%app.SavedData.BestLoopCount    
% for i = [3,7,8,10]%app.SavedData.BestLoopCount    
allClusterMeans = [];
for i = 1:app.SavedData.BestLoopCount
    thisLoopIDs = find(app.SavedData.BestLoopAssignments(:,1) == i);
    thisLoopClusters = app.SavedData.BestLoopAssignments(thisLoopIDs,2);
    
    allIDs = find(app.SavedData.BestStateMap(:,1) == i);
    trialIDs = floor((allIDs-1)/trialLength) + 1;
    conditionID = i;
    
    badIDs = find(app.SavedData.BestStateMap(:,1) == i);
    badIDs = setdiff(1:size(finalStream,1), badIDs);
    
    plotStream = finalStream;
    plotStream(badIDs,:) = nan;
    
    if PLOT_DELAY_EMBEDDED
        h = plot3(plotStream(trialIndicies,:)*trialPCABasis(:,1), plotStream(trialIndicies,:)*trialPCABasis(:,2), plotStream(trialIndicies,:)*trialPCABasis(:,3), 'LineWidth', 1, 'Color', lineColors(conditionID,:));
    else
        h = plot3(plotStream(trialIndicies,1), plotStream(trialIndicies,2), plotStream(trialIndicies,3), 'LineWidth', 1, 'Color', lineColors(conditionID,:));
    end
    h.Color(4) = 0.9;
    
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
        
        allClusterMeans(i, thisLoopClusters(j),:) = clusterMeans(j,:);
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
%     [~, startCluster] = min(meanTimes);
    startCluster = 1;
    
    sortedOrder = [startCluster:length(clusterOrder) 1:startCluster-1];
    thisLoopIDs = thisLoopIDs(sortedOrder);
    
%     for j = 1:length(clusterOrder)
%         thisIndices = find(app.SavedData.BestLoopAssignments(:,1) == i & app.SavedData.BestLoopAssignments(:,2) == clusterOrder(sortedOrder(j)));
%     end
    
    meanTimes = meanTimes(sortedOrder);
    
    badTimes = find(meanTimes <= min(times) | meanTimes >= max(times));
    goodTimes = setdiff(sortedOrder, badTimes);
    
%     thisTrace = app.SavedData.BestEmission(thisLoopIDs(goodTimes), :);
    clusterSTDs = clusterSTDs(thisLoopIDs(goodTimes),:);
    thisTrace = clusterMeans(thisLoopIDs(goodTimes),:);
    
    if length(goodTimes) <= 1
        continue;
    end
    
    thisTrace = filterData(thisTrace', 0.5, [], 1, 0)';
    
    if PLOT_DELAY_EMBEDDED
%         plot3(thisTrace*trialPCABasis(:,1), thisTrace*trialPCABasis(:,2), thisTrace*trialPCABasis(:,3), 'LineWidth', 2, 'Color', colors(i,:))
    
        tubeMean = thisTrace*trialPCABasis(:,1:3);
        projectedSTDs = clusterSTDs*trialPCABasis(:,1:3);
    else
%         plot(thisTrace(:,1), thisTrace(:,2), 'LineWidth', 2, 'Color', colors(i,:))
        
        tubeMean = [thisTrace(:,1:3) zeros(size(thisTrace,1), 1)];
        projectedSTDs = [clusterSTDs(:,1:3) zeros(size(thisTrace,1), 1)];
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
    
    orthogonalSTDs = filterData(orthogonalSTDs, 0.5, [], 1, 2);
    
%     tubeMean = [tubeMean; tubeMean(1,:)];
%     orthogonalSTDs = [orthogonalSTDs; orthogonalSTDs(1,:)];
    
    [X,Y,Z,V] = tubeplot(tubeMean(:,1), tubeMean(:,2), tubeMean(:,3), orthogonalSTDs, ones(1,size(tubeMean,1)),10,[0 0 1]);

    t = surf(X,Y,Z,V, 'FaceColor', lineColors(conditionID,:));
%     if conditionID > 3
%         t.EdgeColor = [0 0 0];
%         t.EdgeAlpha = 0.7;
%     else
%         t.EdgeColor = lineColors(conditionID,:);
%         t.EdgeAlpha = 0.5;
%     end
    t.EdgeColor = lineColors(conditionID,:);
    t.EdgeAlpha = 0.5;
    t.FaceAlpha = 0.2;
end

%% Reconstruction of data

clusterMeans = [];
clusterSTDs = [];
clusterMap = [];
reconstructedData = zeros(size(finalStream));
for i = 1:size(app.SavedData.BestLoopAssignments, 1)
    thisLoopID = app.SavedData.BestLoopAssignments(i,1);
    thisPhase = app.SavedData.BestLoopAssignments(i,2);
    
    thisIndices = find(app.SavedData.BestStateMap(:,1) == thisLoopID & app.SavedData.BestStateMap(:,2) == thisPhase);

    clusterMap(thisIndices) = i;
    
    clusterMeans(i,:) = nanmean(finalStream(thisIndices, :), 1);
    clusterSTDs(i,:) = nanstd(finalStream(thisIndices, :), [], 1);
    
    reconstructedData(thisIndices,:) = repmat(clusterMeans(i,:), [length(thisIndices), 1]);
end
    
    
names = {'x', 'y', 'z'};

figureHandle = figure(2);
figureHandle.Renderer='Painters';
clf;
h = [];
for i = 1:3
    h(i) = subplot(3,1,i);
    hold on;
    plot(finalStream(:,i));
    plot(reconstructedData(:,i));
    ylabel(names{i});
end
xlabel('time');
linkaxes(h,'x');
xlim([1100 1450]);

corrStream = finalStream;
corrStream(isnan(corrStream)) = 0;

reconstructionR2s = [];
for i = 1:size(corrStream,2)
    reconstructionR2s(i) = corr(corrStream(:,i), reconstructedData(:,i));
end
mean(reconstructionR2s)


