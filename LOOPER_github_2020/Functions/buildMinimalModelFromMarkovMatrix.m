

if ~exist('putativeLoopCounts', 'var') || isempty(putativeLoopCounts)
    putativeLoopCounts = [5, 4, 3, 2, 1];
end

if ~exist('totalClusters', 'var') || isempty(totalClusters)
    totalClusters = 100;
end

if ~exist('selectingLoops', 'var') || isempty(selectingLoops)
    selectingLoops = 1;
end

%%

% figure(7);
% clf;
% hold on;
% h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
% h.Color(4) = 0.2;
%
% colors = jet(size(modelMatrix,1));
%
% for j = 1:size(modelMatrix,1)
%     thisID = j;
%
%     points(1,:) = clusterMeansPCA(thisID,:);
%     scatter3(points(:,1), points(:,2), points(:,3), 'o', 'LineWidth', 2,'MarkerEdgeColor',colors(j,:));
% end

%%
%
% cursorMode = datacursormode(gcf);
% datatips = cursorMode.getCursorInfo();
% find(ismember(clusterMeansPCA, datatips(1).Position, 'rows'))

%%

%Add terminal state

modelMatrix = reducedMatrix;
terminalStateID = size(modelMatrix,1) + 1;

if shouldUseTerminalState
    modelMatrix(end+1, end+1) = 1;
end

transitionToMatrix = zeros(size(modelMatrix));
transitionFromMatrix = zeros(size(modelMatrix));
    
transformedTrialEnds = trialSwitchTimes;% - 2*(1:length(trialSwitchTimes));
toTerminalTimes = [transformedTrialEnds, length(clusterIDs)];
fromTerminalTimes = [1, transformedTrialEnds + 1];
fromTerminalTimes(end) = [];

    
for clusterID = 1:size(transitionFromMatrix,1)
    if shouldUseTerminalState && clusterID == terminalStateID        
        fromIDs = toTerminalTimes;
        toIDs = fromTerminalTimes;
    else
        thisIndices = find(clusterIDs == clusterID);

        trajectoryData = regionprops(clusterIDs == clusterID, 'PixelIdxList');

        trajectoryList = {};
        startIDs = [];
        endIDs = [];
        for i = 1:length(trajectoryData)
            trajectoryList{i} = trajectoryData(i).PixelIdxList;

            startIDs(i) = trajectoryList{i}(1);
            endIDs(i) = trajectoryList{i}(end);
        end

        fromIDs = startIDs - 1;
        toIDs = endIDs + 1;
        if shouldUseTerminalState
            fromIDs(ismember(startIDs, fromTerminalTimes)) = 0;
            toIDs(ismember(endIDs, toTerminalTimes)) = 0;
        else
            fromIDs(ismember(startIDs, fromTerminalTimes)) = [];
            toIDs(ismember(endIDs, toTerminalTimes)) = [];
        end
    end

    fromClusters = [];
    toClusters = [];
    if shouldUseTerminalState
        fromClusters = ones(1, sum(fromIDs == 0)) * terminalStateID;
        toClusters = ones(1, sum(toIDs == 0)) * terminalStateID;
    end
    
    fromIDs(fromIDs == 0) = [];
    toIDs(toIDs == 0) = [];
    
    if max(fromIDs) > length(clusterIDs)
        test = 1;
    end
    fromClusters = [fromClusters'; clusterIDs(fromIDs)];
    toClusters = [toClusters'; clusterIDs(toIDs)];
    
    for i = unique(toClusters)'
        transitionToMatrix(clusterID, i) = sum(toClusters == i);
    end
    %     transitionToMatrix(clusterID,:) = transitionToMatrix(clusterID,:) / sum(transitionToMatrix(clusterID,:));
    
    for i = unique(fromClusters)'
        transitionFromMatrix(clusterID, i) = sum(fromClusters == i);
    end
    %     transitionFromMatrix(clusterID,:) = transitionFromMatrix(clusterID,:) / sum(transitionFromMatrix(clusterID,:));
    
end

%%

MIN_VALUE = 1/10000;
DISPLAY = 1;

graphTree = (transitionToMatrix);
graphTree(graphTree < MIN_VALUE) = MIN_VALUE;
graphTree(find(eye(size(graphTree)))) = MIN_VALUE;
graphTree = 1 ./ graphTree;

% if DISPLAY
%     figure(5);
%     clf;
%     hold on;
%     h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
%     h.Color(4) = 0.2;
% end

allLoops  = {};
loopCorrelationMatrix = zeros(size(transitionToMatrix,1));
loopValues = [];
loopCounts = zeros(1, size(transitionToMatrix,1));
for i = 1:size(transitionToMatrix,1)
    loopStart = i;
    possibleEnds = find(transitionFromMatrix(loopStart,:));
    bestLength = 99999999;
    bestPath = [];
    for j = 1:length(possibleEnds)
        if possibleEnds(j) == loopStart
            continue;
        end
        
        if shouldUseTerminalState
            [pathLength, path, pred] = graphshortestpath(sparse(graphTree), loopStart, possibleEnds(j));
            
            [pathLength1, path1, pred1] = graphshortestpath(sparse(graphTree), loopStart, terminalStateID);
            [pathLength2, path2, pred2] = graphshortestpath(sparse(graphTree), terminalStateID, possibleEnds(j));
            terminalPathLength = pathLength1 + pathLength2;
            if length(path1) > 1
               terminalPathLength = terminalPathLength - (graphTree(path1(end-1), terminalStateID));
            end
            if length(path2) > 1
               terminalPathLength = terminalPathLength - (graphTree(terminalStateID, path2(2)));
            end
            terminalPath = [path1(1:end-1) path2];
            
            if terminalPathLength < pathLength
                pathLength = terminalPathLength;
                path = terminalPath;
            end
        else
            [pathLength, path, pred] = graphshortestpath(sparse(graphTree), loopStart, possibleEnds(j));
        end
        
        totalLength = pathLength + graphTree(possibleEnds(j), loopStart);
        totalLength = totalLength ./ sqrt(length(path));
        
        DISPLAY=1;
        if DISPLAY && i == 1
            figure(7);
            clf;
            hold on;
            h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
            h.Color(4) = 0.2;
            
            for j = 1:size(clusterMeansPCA,1)
                points(1,:) = clusterMeansPCA(j,:);
                points(2,:) = clusterMeansPCA(j,:) + 2;
                
                text(points(2,1), points(2,2), points(2,3), num2str(j), 'horizontal', 'center')
                
                
                h = plot3(points(:,1), points(:,2), points(:,3), 'k-', 'LineWidth', 0.5);
                h.Color(4) = 0.5;
            end
            
            scatter3(clusterMeansPCA(:,1), clusterMeansPCA(:,2), clusterMeansPCA(:,3), 'b', 'filled')
            axis off
            path;
            pathNoise = normrnd(0, 0.2, [length(path), size(clusterMeansPCA(1,:),2)]);
            pathNoise(end,:) = pathNoise(1,:);
            points = [];
            for j = 1:length(path)-1
                thisID = path(j);

                if thisID > size(clusterMeansPCA,1) || path(j+1)  > size(clusterMeansPCA,1)
                    continue;
                end

                points(1,:) = clusterMeansPCA(thisID,:) + pathNoise(j,:);
                points(2,:) = clusterMeansPCA(path(j+1),:) + pathNoise(j+1,:);
                plot3(points(:,1), points(:,2), points(:,3), '-ro', 'LineWidth', 2);
            end
            
            title(num2str(i))
            
            pause(1);
        end
        
        %         if length(bestPath) < 3 || (totalLength < bestLength && length(path) > 2)
        if totalLength < bestLength
            bestLength = totalLength;
            bestPath = path;
        end
    end
    
    allLoops{loopStart} = bestPath;
    loopCorrelationMatrix(loopStart,bestPath) = 1;
    loopValues(loopStart) = bestLength;
    
    %     if DISPLAY && i == uniqueLoops(10)%i >= 26 && i <= 26
    %         colors = lines(size(transitionToMatrix,1));
    %         drawPath = [bestPath bestPath(1)];
    %         for j = 1:length(drawPath)-1
    %             thisID = drawPath(j);
    %
    %             points(1,:) = clusterMeansPCA(thisID,:);
    %             points(2,:) = clusterMeansPCA(drawPath(j+1),:);
    %             plot3(points(:,1), points(:,2), points(:,3), '-o', 'LineWidth', 2,'Color',colors(loopStart,:));
    %         end
    %     end
    
end


loopBandwidth = [];
loopWeight = [];
for i = 1:size(transitionToMatrix,1)
    loopIDs = allLoops{i};
    
    transitionBandwidth = [];
    transitionWeight = [];
    for j = 1:length(loopIDs)
        nextID = j + 1;
        if nextID > length(loopIDs)
            nextID = 1;
        end
        
        transitionBandwidth(j) = transitionToMatrix(loopIDs(j), loopIDs(nextID));
        transitionWeight(j) = modelMatrix(loopIDs(j), loopIDs(nextID));
    end
    
    loopBandwidth(i) = mean(transitionBandwidth);
    loopWeight(i) = sum(transitionWeight);
end
loopBandwidth(loopBandwidth == 0) = nan;

%%

% [steadyState, ~] = eigs(asymmetricProbabilities', 1);
% fluxMatrix = diag(steadyState) * asymmetricProbabilities;
% symmetricMatrix = (fluxMatrix+fluxMatrix.')/2;
% antisymmetricMatrix = (fluxMatrix-fluxMatrix.')/2;
% symmetricMatrix = symmetricMatrix ./ steadyState;
% antisymmetricMatrix = antisymmetricMatrix ./ steadyState;
% 
% stateTransitions = [];
% for i = 1:size(modelMatrix,1)
%     thisStates = find(clusterIDs == i);
%     
%     stateTransitions(i,:) = sum(symmetricMatrix(thisStates,:));
%     stateTransitions(i,:) = stateTransitions(i,:) ./ norm(stateTransitions(i,:));
% end
% 
% stateSimilarities = stateTransitions * stateTransitions';

%%

loopOverlap = loopCorrelationMatrix * loopCorrelationMatrix';
loopLengths = sum(loopCorrelationMatrix,2);
loopOverlap = loopOverlap ./ loopLengths;

mergedLoopIDs = 1:size(modelMatrix,1);
removedLoops = [];
for i = 1:size(loopCorrelationMatrix,1)
    overlapLoops = find(loopOverlap(i,:) == 1);
    overlapLoops = sort(setdiff(overlapLoops, [removedLoops, i]), 'descend');
    [~, longestLoop] = max(loopLengths(overlapLoops));
    
    if length(longestLoop) > 0
        mergeLoop = overlapLoops(longestLoop);
        
        mergedLoopIDs(i) = mergeLoop;
        
        if mergeLoop ~= i
            removedLoops = [removedLoops, i];
        end
    end
end

mergedLoopIDs = 1:size(loopCorrelationMatrix,1);

uniqueLoops = unique(mergedLoopIDs);
uniqueCorrelations = loopCorrelationMatrix(uniqueLoops, :);
uniqueBandwidths = loopBandwidth(uniqueLoops);

otherStateSimilarities = stateSimilarities .* (1 - eye(size(stateSimilarities)));

loopSimilarities = [];
for i = 1:length(uniqueLoops)
    for j = 1:length(uniqueLoops)
        mismatches = find(uniqueCorrelations(i,:) == 1 & uniqueCorrelations(j,:) == 0);
        mismatchesOther = find(uniqueCorrelations(i,:) == 0 & uniqueCorrelations(j,:) == 1);
        %         mismatchesOther = find(uniqueCorrelations(j,:) == 1);
        
        mismatches(mismatches > size(otherStateSimilarities,1)) = [];
        mismatchesOther(mismatchesOther > size(otherStateSimilarities,1)) = [];
        
        if length(mismatches) > 0 && length(mismatchesOther) > 0
            mismatchSimilarities = [];
            for k = mismatches
                mismatchSimilarities(end+1) = max(otherStateSimilarities(k,mismatchesOther));
            end
            
            loopSimilarities(i,j) = (prod(mismatchSimilarities));
        elseif length(mismatches) == 0
            loopSimilarities(i,j) = 1;
        else
            loopSimilarities(i,j) = 0;
        end
    end
end

weightedLoopSimilarities = loopSimilarities ./ uniqueBandwidths;%' * uniqueBandwidths);
weightedLoopSimilarities = log(1 ./ max(weightedLoopSimilarities, weightedLoopSimilarities'));
% weightedLoopSimilarities = 1 ./ m(weightedLoopSimilarities * weightedLoopSimilarities');
weightedLoopSimilarities = weightedLoopSimilarities .* (1 - eye(size(weightedLoopSimilarities)));

weightedLoopSimilarities(isnan(weightedLoopSimilarities)) = 0;

clusterSimilarities = squareform(weightedLoopSimilarities);

loopClustering = linkage(clusterSimilarities, 'complete');


%% Check loop clusters

TOTAL_LOOPS = totalClusters;

allLikelihoods = [];
% averageVariances = [];
allModels = {};
allEmissions = {};
allLoopAssignments = {};
allStateMaps = {};


if ~selectingLoops
    waitHandle = parfor_progressbar(length(putativeLoopCounts), 'Calculating model accuracy');
end
for maxClusterIndex = 1:length(putativeLoopCounts)
    maxClusters = putativeLoopCounts(maxClusterIndex);
    
    loopClusterIDs = cluster(loopClustering,'maxclust',maxClusters);
    
    figure(6);
    clf
    cutoff = median([loopClustering(end-maxClusters+1,3) loopClustering(end-maxClusters+2,3)]);
    dendrogram(loopClustering,'ColorThreshold',cutoff)
    
    if ~exist('app')
        app = [];
    end
    UIData = drawLoops(app, selectingLoops, false, randi(1), maxClusters, loopClusterIDs, finalDynamicsStream, pcaBasis, allLoops, uniqueLoops, clusterMeansPCA, loopBandwidth, mergedLoopIDs, TOTAL_LOOPS, clusterIDs, loopCorrelationMatrix, stateSimilarities, loopWeight, loopLengths, clusterMeans, countClusters);
    
    if selectingLoops        
        
        
                     test=1;
    else
        buildMinimalModelFromLoops
        
        %%
        
        validationModel = minimalModel;
        validationEmission = statePositions;
        
        validateModel;
        
        allLikelihoods(end+1) = scoreMean;
        
        allModels{end+1} = minimalModel;
        allEmissions{end+1} = statePositions;
        allLoopAssignments{end+1} = modelLoopAssignments;
        allStateMaps{end+1} = stateMap;
        
        if maxClusters == 3
            tes = 1;
        end
        
        waitHandle.iterate(1);
    end
end
if ~selectingLoops
    close(waitHandle);
end

%%

if ~selectingLoops
    
    allLikelihoods
    
    [~, bestLoopID] = min(allLikelihoods);
    % bestLoopID = 3;
    
    bestLoopCount = putativeLoopCounts(bestLoopID);
    
    bestModel = allModels{bestLoopID};
    bestEmission = allEmissions{bestLoopID};
    bestLoopAssignments = allLoopAssignments{bestLoopID};
    bestStateMap = allStateMaps{bestLoopID};
    
    bestStateCount = size(bestModel,1);
    
    % totalParameters = putativeLoopCounts / 2 * log(numDataPoints / (2*pi));
    %
    %
    % BICs = totalLikelihoods' + totalParameters;
    
    
end


%callback function for the pushbutton  (save it in its own *.m file if needed)           
function buttonSelected(varargin)
    UIData = varargin{3};  % Get the structure.
    
    if UIData.cancel == varargin{1}
        buttonID = 0;
    elseif UIData.done == varargin{1}
        UIData.app.updateLoops(UIData);
        
        return;
    else
        buttonID = find(UIData.buttons == varargin{1});

        if UIData.currentSelection > 0
            UIData.maxClusters = UIData.maxClusters - 1;

            toID = min(UIData.currentSelection, buttonID);
            fromID = max(UIData.currentSelection, buttonID);

            UIData.loopClusterIDs(UIData.loopClusterIDs == fromID) = toID;
            UIData.loopClusterIDs(UIData.loopClusterIDs > fromID) = UIData.loopClusterIDs(UIData.loopClusterIDs > fromID) - 1;
            test = 1;

            buttonID = 0;
        end

    %     delete(UIData.buttons(buttonID))
    end
    
    drawLoops(UIData.app, true, buttonID, UIData.seed, UIData.maxClusters, UIData.loopClusterIDs, UIData.finalDynamicsStream, UIData.pcaBasis, UIData.allLoops, UIData.uniqueLoops, UIData.clusterMeansPCA, UIData.loopBandwidth, UIData.mergedLoopIDs, UIData.TOTAL_LOOPS, UIData.clusterIDs, UIData.loopCorrelationMatrix, UIData.stateSimilarities, UIData.loopWeight, UIData.loopLengths, UIData.clusterMeans, UIData.countClusters);
end

function [UIData] = drawLoops(app, makeButtons, isSelected, seed, maxClusters, loopClusterIDs, finalDynamicsStream, pcaBasis, allLoops, uniqueLoops, clusterMeansPCA, loopBandwidth, mergedLoopIDs, TOTAL_LOOPS, clusterIDs, loopCorrelationMatrix, stateSimilarities, loopWeight, loopLengths, clusterMeans, countClusters)
    UIData = struct;
    
    if makeButtons
        UIData.app = app;
        
        if isSelected == 0
            rng shuffle
            UIData.seed = randi(1);
        else 
            UIData.seed = seed;
        end

        rng(UIData.seed);

        UIData.currentSelection = isSelected;
        UIData.maxClusters = maxClusters;
        UIData.loopClusterIDs = loopClusterIDs;
        UIData.finalDynamicsStream = finalDynamicsStream;
        UIData.pcaBasis = pcaBasis;
        UIData.allLoops = allLoops;
        UIData.uniqueLoops = uniqueLoops;
        UIData.clusterMeans = clusterMeans;
        UIData.clusterMeansPCA = clusterMeansPCA;
        UIData.loopBandwidth = loopBandwidth;
        UIData.mergedLoopIDs = mergedLoopIDs;
        UIData.TOTAL_LOOPS = TOTAL_LOOPS;
        UIData.clusterIDs = clusterIDs;
        UIData.loopCorrelationMatrix = loopCorrelationMatrix;
        UIData.stateSimilarities = stateSimilarities;
        UIData.loopWeight = loopWeight;
        UIData.loopLengths = loopLengths;
        UIData.countClusters = countClusters;
        
    end

    xSize = ceil(sqrt(maxClusters));
    ySize = ceil(maxClusters / xSize);
    UIData.f = figure(5);
    clf;
    if isSelected == 0
        drawnow;
    end
    for l = 1:maxClusters
        plotLoops = uniqueLoops(find(loopClusterIDs == l));
        
        subh = subplot(xSize, ySize, l);
        hold on;
        h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
        h.Color(4) = 0.2;
        
        axis off
        
        colors = parula(max(round(loopBandwidth)));
        loopIndex = 1;
        for i = 1:length(plotLoops)
            allClusteredLoops = find(mergedLoopIDs == plotLoops(i));
            for k = 1:length(allClusteredLoops)
                loopID = allClusteredLoops(k);
                
                if isempty(allLoops{loopID})
                    continue;
                end
                
                drawPath = [allLoops{loopID} allLoops{loopID}(1)];
                pathNoise = normrnd(0, 0.2, [length(drawPath), size(clusterMeansPCA(1,:),2)]);
                pathNoise(end,:) = pathNoise(1,:);
                points = [];
                for j = 1:length(drawPath)-1
                    thisID = drawPath(j);
                    
                    if thisID > size(clusterMeansPCA,1) || drawPath(j+1)  > size(clusterMeansPCA,1)
                        continue;
                    end
                    
                    points(1,:) = clusterMeansPCA(thisID,:) + pathNoise(j,:);
                    points(2,:) = clusterMeansPCA(drawPath(j+1),:) + pathNoise(j+1,:);
                    plot3(points(:,1), points(:,2), points(:,3), '-o', 'LineWidth', 2,'Color',colors(round(loopBandwidth(loopID)),:));
                end
                
                loopIndex = loopIndex + 1;
            end
        end
        
        if makeButtons
            text = 'Merge this loop';
            enabled = 'off';
            if isSelected ~= l
                enabled = 'on';
            end
            if isSelected > 0 && isSelected ~= l
                text = 'Into this loop';
            end
            
            UIData.buttons(l) = uicontrol('style','push',...
                     'units','pix',...
                     'position', [subh.Position(1) * UIData.f.Position(3) - 20, subh.Position(2) * UIData.f.Position(4) - 20, 150, 25],...
                     'fontsize',14,...
                     'string', text,...
                 'enable', enabled);
        end
    end
    
    if makeButtons
        enabled = 'off';
        if isSelected > 0
            enabled = 'on';
        end
        UIData.cancel = uicontrol('style','push',...
                     'units','pix',...
                     'position', [160, 10, 150, 25],...
                     'fontsize',14,...
                     'string', 'Cancel',...
                 'enable', enabled);
             
        enabled = 'on';
        if isSelected > 0
            enabled = 'off';
        end
         UIData.done = uicontrol('style','push',...
                     'units','pix',...
                     'position', [10, 10, 150, 25],...
                     'fontsize',14,...
                     'string', 'Finish',...
                 'enable', enabled);
        
        for i = 1:length(UIData.buttons)
            set(UIData.buttons(i), 'callback', {@buttonSelected,UIData});
        end
        
        set(UIData.cancel, 'callback', {@buttonSelected,UIData});
        set(UIData.done, 'callback', {@buttonSelected,UIData});
    end
    
end

