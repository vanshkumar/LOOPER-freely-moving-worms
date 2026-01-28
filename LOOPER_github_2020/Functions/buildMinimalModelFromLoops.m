
%% Align loops

BINS_PER_LOOP = ceil(TOTAL_LOOPS / maxClusters);
FIG_BASE = 400;

terminalStateID = max(clusterIDs) + 1;

averageLoops = [];
averageErrors = [];
totalPoints = [];
stateMap = zeros(size(clusterIDs,1),2);
for i = 1:maxClusters
    clusterLoops = uniqueLoops(find(loopClusterIDs == i));
    
    if ismember(terminalStateID, clusterLoops)
        hasTerminalState = 1;
    else
        hasTerminalState = 0;
    end

    if length(clusterLoops) > 1
        [~, bestStartClusterID] = max(sum(loopCorrelationMatrix(clusterLoops, clusterLoops)));
        bestStartClusterID = clusterLoops(bestStartClusterID);
    else
        bestStartClusterID = clusterLoops;
    end
    
    [~, bestLoopID] = max(loopBandwidth(clusterLoops));
    bestLoopID = clusterLoops(bestLoopID);

    loopPhases = {};
    loops = {};
    clusterPhases = nan(length(clusterLoops), length(loopClusterIDs));
    
    loopLength = length(allLoops{bestLoopID});
    path = allLoops{bestLoopID};
    [~, startPosition] = max(stateSimilarities(bestStartClusterID, path));

    loopPhases{1} = ((1:loopLength) - 1)/(loopLength+1)*2*pi;
    loops{1} = circshift(allLoops{bestLoopID}, -(startPosition-1));
    clusterPhases(1, loops{end}) = loopPhases{end};

    for j = 1:length(clusterLoops)
        loopID = clusterLoops(j);

        if loopID == bestLoopID
            continue;
        end

        loopLength = length(allLoops{loopID});
        path = allLoops{loopID};
        [~, startPosition] = max(stateSimilarities(bestStartClusterID, path));

        loopPhases{end+1} = ((1:loopLength) - 1)/(loopLength+1)*2*pi;
        loops{end+1} = circshift(allLoops{loopID}, -(startPosition-1));

        clusterPhases(j, loops{end}) = loopPhases{end};
    end

    thisLoopWeights = loopWeight(clusterLoops);
    thisLoopWeights = thisLoopWeights / sum(thisLoopWeights);
    ANGLE_SIGMA = 2*pi/(sum(loopLengths(clusterLoops).*thisLoopWeights')*2);
    if isnan(ANGLE_SIGMA)
        ANGLE_SIGMA = 1;
    end
    
    if hasTerminalState
        terminalStatePhase = nanmean(clusterPhases(:,terminalStateID));
        clusterPhases = clusterPhases - terminalStatePhase;
        
        clusterPhases(:,1:end-1) = wrapTo2Pi(clusterPhases(:,1:end-1));
    end
    
    
    
    loopPositions = zeros(BINS_PER_LOOP, size(clusterMeans,2));
    loopPhase = [];
    for j = 1:BINS_PER_LOOP
        testPhase = (j-1) / (BINS_PER_LOOP+1) * 2*pi;
        
        if hasTerminalState
            phaseWeights = exp(-(abs(testPhase - clusterPhases)).^2 / (2*ANGLE_SIGMA^2));
        else
            phaseWeights = exp(-(angleDiff(testPhase, clusterPhases)).^2 / (2*ANGLE_SIGMA^2));
        end

        phaseWeights = phaseWeights .* loopWeight(clusterLoops)';

        if shouldUseTerminalState
            phaseWeights = nansum(phaseWeights) .* sqrt([countClusters 0]);
        else
            phaseWeights = nansum(phaseWeights) .* sqrt(countClusters);
        end
        
        if nansum(phaseWeights) == 0
            phaseWeights(clusterLoops) = 1;
        else
            phaseWeights = phaseWeights / nansum(phaseWeights);
        end

        if shouldUseTerminalState
            loopPositions(j,:) = phaseWeights(1:end-1) * clusterMeans;
        else
            loopPositions(j,:) = phaseWeights * clusterMeans;
        end
        loopPhase(j) = testPhase;
    end

    [loopPCABasis, ~, ~, ~, explained] = pca(loopPositions);
    requiredDimensions = find(cumsum(explained) > 99, 1, 'first');
    loopPCABasis = loopPCABasis(:, 1:requiredDimensions);

    allClusters = [];
    for j = 1:length(clusterLoops)
        allClusters = unique([allClusters, allLoops{clusterLoops(j)}]);
    end

    thisPoints = [];
    for j = 1:length(allClusters)
        thisPoints = unique([thisPoints, find(clusterIDs == allClusters(j))']);
    end

    
    figure(FIG_BASE + 5);
    clf;
    hold on;
    if size(pcaBasis,2) > 2
    	h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
    else
        h = plot(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), 'LineWidth', 0.5);
    end
    h.Color(4) = 0.2;
    plotData = loopPositions;
    plotData = plotData + normrnd(0, 0.2, size(plotData));
    if size(pcaBasis,2) > 2
        plot3(plotData*pcaBasis(:,1), plotData*pcaBasis(:,2), plotData*pcaBasis(:,3), '-', 'LineWidth', 1);
    else
        plot(plotData*pcaBasis(:,1), plotData*pcaBasis(:,2), '-', 'LineWidth', 1);
    end

    pointDistances = pdist2(finalDynamicsStream(thisPoints,:)*loopPCABasis, loopPositions*loopPCABasis);

    title('Loop positions over dynamics (PCA space)');
    if size(pcaBasis,2) <= 2
        xlabel('PC1');
        ylabel('PC2');
    end

    [~, pointLoopPosition] = min(pointDistances, [], 2);

    DISPLAY = 0;
    if DISPLAY
        colors = hsv(BINS_PER_LOOP);
        testLoopColors = hsv(size(pointDistances,2));
    
        figure(FIG_BASE + 4);
        clf;
        hold on
        if size(pcaBasis,2) > 2
            h = plot3(finalDynamicsStream*loopPCABasis(:,1), finalDynamicsStream*loopPCABasis(:,2), finalDynamicsStream*loopPCABasis(:,3), 'LineWidth', 0.5);
        else
            h = plot3(finalDynamicsStream*loopPCABasis(:,1), finalDynamicsStream*loopPCABasis(:,2), 'LineWidth', 0.5);
        end
        h.Color(4) = 0.2;
        
        plotData = [loopPositions; loopPositions(1,:)];
        if size(pcaBasis,2) > 2
            plot3(plotData*loopPCABasis(:,1), plotData*loopPCABasis(:,2), plotData*loopPCABasis(:,3), 'k-', 'LineWidth', 1);
        else
            plot(plotData*loopPCABasis(:,1), plotData*loopPCABasis(:,2), 'k-', 'LineWidth', 1);
        end
        plotData = loopPositions;
        if size(pcaBasis,2) > 2
            scatter3(plotData*loopPCABasis(:,1), plotData*loopPCABasis(:,2), plotData*loopPCABasis(:,3), 72, testLoopColors);
        else
            scatter(plotData*loopPCABasis(:,1), plotData*loopPCABasis(:,2), 72, testLoopColors);
        end
    end
    
    finalLoopPositions = [];
    loopError = [];
    loopCounts = [];
    for j = 1:BINS_PER_LOOP
        thisIndices = thisPoints(find(pointLoopPosition == j));
        thisWeights = loopWeight(clusterIDs(thisIndices));
        thisWeights = thisWeights / sum(thisWeights);
        thisWeights(isnan(thisWeights)) = 0;

        stateMap(thisIndices,:) = repmat([i, j], length(thisIndices), 1);

        finalLoopPositions(j,:) = mean(finalDynamicsStream(thisIndices,:), 1);

        loopError(j) = var(pdist2(single(finalDynamicsStream(thisIndices,:)*loopPCABasis / norm(loopPCABasis)), single(finalLoopPositions(j,:)*loopPCABasis / norm(loopPCABasis))), thisWeights);

        loopCounts(j) = length(thisIndices);
        
        if DISPLAY
            hold on;
            if size(pcaBasis,2) > 2
                plot3(finalDynamicsStream(thisIndices,:)*loopPCABasis(:,1), finalDynamicsStream(thisIndices,:)*loopPCABasis(:,2), finalDynamicsStream(thisIndices,:)*loopPCABasis(:,3), 'o', 'LineWidth', 1, 'Color', colors(j,:));
            else
                 plot(finalDynamicsStream(thisIndices,:)*loopPCABasis(:,1), finalDynamicsStream(thisIndices,:)*loopPCABasis(:,2), 'o', 'LineWidth', 1, 'Color', colors(j,:));
            end
        end
    end
    
    if hasTerminalState
        finalLoopPositions = filterData(finalLoopPositions', 0.5, 'gaussian', 1, 0)';
    else
        finalLoopPositions = filterData(finalLoopPositions', 0.5, 'gaussian', 1, 2)';
    end
    
    averageLoops(i,:,:) = finalLoopPositions;
    averageErrors(i) = nansum(loopError .* loopCounts / sum(loopCounts));
    totalPoints(i) = length(thisPoints);
end

totalError = sum(averageErrors .* totalPoints / sum(totalPoints));


[uniqueStates, ~, stateIDs] = unique(stateMap, 'row');

modelLoopAssignments = uniqueStates;

statePositions = [];
for i = 1:size(uniqueStates,1)    
    statePositions(i,:) = averageLoops(uniqueStates(i,1), uniqueStates(i,2),:);
end

minimalModel = zeros(size(uniqueStates,1),size(uniqueStates,1));
for i = 1:length(stateIDs)-1
    thisState = stateIDs(i);
    nextState = stateIDs(i+1);

    minimalModel(thisState, nextState) = minimalModel(thisState, nextState) + 1;
end
minimalModel = minimalModel ./ sum(minimalModel,2);

DEBUG = 0;
if DEBUG == 1
    figure(6);
    clf;
    hold on;
    if size(pcaBasis,2) > 2
        h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
    else
        h = plot(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), 'LineWidth', 0.5);
    end
    h.Color(4) = 0.2;
    plotData = statePositions(stateIDs,:);
    plotData = plotData + normrnd(0, 0.2, size(plotData));
    if size(pcaBasis,2) > 2
        plot3(plotData*pcaBasis(:,1), plotData*pcaBasis(:,2), plotData*pcaBasis(:,3), '-', 'LineWidth', 1);
    else
        plot(plotData*pcaBasis(:,1), plotData*pcaBasis(:,2), '-', 'LineWidth', 1);
    end
end
