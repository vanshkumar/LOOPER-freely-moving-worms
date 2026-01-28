
%%

MAX_COMMUNITY_ATTEMPTS = 10;
USE_PHASE_VELOCITY = 2*pi/100;
SMOOTH_PHASE = 2*pi/20;
END_SEGMENTS = 0;%round(SMOOTH_PHASE/USE_PHASE_VELOCITY);

[eigenVectors, eigenValues] = eigs(asymmetricProbabilities, 10);
eigenValues = diag(eigenValues);

[~, bestComplexEigenVector] = max((abs(imag(eigenValues)) > 0) .* abs(eigenValues));

globalPhase = angle(eigenVectors(:,bestComplexEigenVector));

trajectoryMeans = {};
trajectorySTDs = {};
trajectoryWeights = {};
trajectoryClusterMap = [];
trajectoryQuality = {};
trajectoryCertainty = {};
trajectoryCounts = [];
trajectoryIndices = {};
trajectoryPhases = [];
dataTrajectoryMap = [];

if ~exist('DISPLAY')
    DISPLAY = 1;
end
if ~exist('DISPLAY_INCREMENT')
    DISPLAY_INCREMENT = 0;
end

if ~DISPLAY_INCREMENT
    waitHandle = parfor_progressbar(size(reducedMatrix,1), 'Building local trajectories');
end
for clusterID = 1:size(reducedMatrix,1)
    % for clusterID = 48
    % for clusterID = 45
    % for clusterID = 27
    %%
    %     if ~validClusters(clusterID)
    %     	continue
    %     end
    
    thisIndices = find(clusterIDs == clusterID);
    
    %     figure(1);
    %     clf;
    %     hold on;
    %     h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
    %     h.Color(4) = 0.2;
    %     scatter3(finalDynamicsStream(thisIndices,:)*pcaBasis(:,1), finalDynamicsStream(thisIndices,:)*pcaBasis(:,2), finalDynamicsStream(thisIndices,:)*pcaBasis(:,3), 32, phase(i)*ones(size(finalDynamicsStream(thisIndices,3))));
    
    %%
    
    MIN_CLUSTER_COUNT = 2;
    TAIL_LENGTH = 20;
    
    trajectoryList = {};
    inGroupTrajectoryList = {};
    trajectoryIDs = {};
    inGroupTrajectoryIDs = {};
    
    startIDs = [];
    endIDs = [];
    startClusters = [];
    endClusters = [];
    groupIndices = [];
    currentID = 0;
    trajectoryData = regionprops(clusterIDs == clusterID, 'PixelIdxList');
    
    if length(trajectoryData) < MIN_CLUSTER_COUNT
        continue;
    end
    
    for i = 1:length(trajectoryData)
        inGroupTrajectoryList{i} = trajectoryData(i).PixelIdxList;
        
        tempIndices = trajectoryData(i).PixelIdxList + (-TAIL_LENGTH:TAIL_LENGTH);
        tempIndices(tempIndices < 1) = [];
        tempIndices(tempIndices > length(clusterIDs)) = [];
        tempIndices = unique(tempIndices);
        if size(tempIndices,2) > 1
            tempIndices = tempIndices';
        end
        
        %         sum(stateHasNext(trajectoryData(i).PixelIdxList) == 0)
        
        tailIndices = setdiff(tempIndices, inGroupTrajectoryList{i});
        startTail = tailIndices;
        endTail = tailIndices;
        startTail(startTail > trajectoryData(i).PixelIdxList(1)) = [];
        endTail(startTail < trajectoryData(i).PixelIdxList(1)) = [];
        
        startCutoff = find(stateHasNext(startTail) == 0, 1, 'last') + 1;
        endCutoff = find(stateHasNext(endTail) == 0, 1) - 1;
        if isempty(startCutoff)
            startCutoff = 1;
        end
        if isempty(endCutoff)
            endCutoff = length(endTail);
        end
        tempIndices = inGroupTrajectoryList{i};
        if startCutoff < length(startTail)
            tempIndices = [startTail(startCutoff:end); tempIndices];
        end
        if endCutoff > 0
            tempIndices = [tempIndices; endTail(1:endCutoff)];
        end
        %
        %         tempIndices = [startTail; inGroupTrajectoryList{i}; endTail];
        
        trajectoryList{i} = tempIndices;
        
        startIDs(i) = inGroupTrajectoryList{i}(1);
        endIDs(i) = inGroupTrajectoryList{i}(end);
        
        startClusters(i) = clusterIDs(max(1, trajectoryData(i).PixelIdxList(1) - 1));
        endClusters(i) = clusterIDs(min(length(clusterIDs), trajectoryData(i).PixelIdxList(end) + 1));
        
    end
    
    for i = 1:length(trajectoryData)
        groupIndices = unique([groupIndices; trajectoryList{i}]);
    end
    
    for i = 1:length(trajectoryData)
        trajectoryIDs{i} = find(ismember(groupIndices, trajectoryList{i}));
        inGroupTrajectoryIDs{i} = find(ismember(groupIndices, inGroupTrajectoryList{i}));
    end
    
    tailIndices = setdiff(groupIndices, thisIndices);
    inGroupIndices = ~ismember(groupIndices, tailIndices);
    
    reducedProbabilities = asymmetricProbabilities(groupIndices,groupIndices);
    reducedSums = sum(reducedProbabilities,2);
    badIndicies = find(reducedSums == 0);
    reducedProbabilities(badIndicies,:) = reducedProbabilities(badIndicies-1,:);
    reducedProbabilities = reducedProbabilities ./ nansum(reducedProbabilities,2);
    
    %     MAX_PROBABILITY = exp(-1);
    %
    %     maxProbabilities = max(reducedProbabilities,[],2);
    %     badIndiciesRow = maxProbabilities > MAX_PROBABILITY;
    %     maxProbabilities = max(reducedProbabilities,[],1);
    %     badIndiciesColumn = maxProbabilities > MAX_PROBABILITY;
    %     badIndicies = unique([find(badIndiciesRow); find(badIndiciesColumn)']);
    %
    %     reducedProbabilities(badIndicies,:) = [];
    %     reducedProbabilities(:,badIndicies) = [];
    
    reducedProbabilities = reducedProbabilities ./ nansum(reducedProbabilities,2);
    
    %% Get detailed balance decomp
    [steadyState, ~] = eigs(reducedProbabilities', 1);
    fluxMatrix = diag(steadyState) * reducedProbabilities;
    symmetricMatrix = (fluxMatrix+fluxMatrix.')/2;
    antisymmetricMatrix = (fluxMatrix-fluxMatrix.')/2;
    symmetricMatrix = symmetricMatrix ./ steadyState;
    antisymmetricMatrix = antisymmetricMatrix ./ steadyState;
    
    %% Rescale diffusion map
    %         diffusedProbabilities = max(0, expm(symmetricMatrix) - eye(size(symmetricMatrix)) + antisymmetricMatrix);
    diffusedProbabilities = max(0, symmetricMatrix + sign(antisymmetricMatrix).*abs(antisymmetricMatrix).^(1/4));
    %     diffusedProbabilities = antisymmetricMatrix;
    %     diffusedProbabilities = max(0, symmetricMatrix + antisymmetricMatrix);
    % diffusedProbabilities = symmetricMatrix + antisymmetricMatrix;
    diffusedProbabilities = diffusedProbabilities ./ sum(diffusedProbabilities,2);
    
    %%
    
    [eigenVectors, eigenValues] = eigs(diffusedProbabilities, 10);
    eigenValues = diag(eigenValues);
    
    [~, bestComplexEigenVector] = max((abs(imag(eigenValues)) > 0) .* abs(eigenValues));
    
    
    localPhase = angle(eigenVectors(:,bestComplexEigenVector));
    %     localPhase = globalPhase(groupIndices);
    meanPhase = angle(sum(exp(1i*localPhase(inGroupIndices))));
    localPhase = wrapToPi(localPhase - meanPhase);
    
    %     for i = 1:length(inGroupTrajectoryIDs)
    %         [~, midPhaseIndex] = min(abs(localPhase(inGroupTrajectoryIDs{i})));
    %
    %         fullTrajectoryID = find(trajectoryIDs{i} == inGroupTrajectoryIDs{i}(midPhaseIndex));
    %
    %         forwardAngles = unwrap(localPhase(trajectoryIDs{i}(fullTrajectoryID:end)));
    %         backwardAngles = unwrap(localPhase(trajectoryIDs{i}(fullTrajectoryID-1:-1:1)));
    %
    %         trajectoryAngles = [backwardAngles(end:-1:1); forwardAngles];
    %
    %         localPhase(trajectoryIDs{i}) = trajectoryAngles;
    %     end
    
    %     figure(2);
    %     clf;
    %     plot(abs(eigenValues));
    %
    
    %     jumpIndicies = find(stateHasNext == 0);
    
    %     figure(3);
    %     clf;
    %     hold on;
    %     h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
    %     h.Color(4) = 0.2;
    %     scatter3(finalDynamicsStream(jumpIndicies,:)*pcaBasis(:,1), finalDynamicsStream(jumpIndicies,:)*pcaBasis(:,2), finalDynamicsStream(jumpIndicies,:)*pcaBasis(:,3), 32, 'kx');
    % %     scatter3(finalDynamicsStream(groupIndices,:)*pcaBasis(:,1), finalDynamicsStream(groupIndices,:)*pcaBasis(:,2), finalDynamicsStream(groupIndices,:)*pcaBasis(:,3), 32, localPhase);
    % %     scatter3(finalDynamicsStream(tailIndices,:)*pcaBasis(:,1), finalDynamicsStream(tailIndices,:)*pcaBasis(:,2), finalDynamicsStream(tailIndices,:)*pcaBasis(:,3), 70, 'kx', 'LineWidth', 0.5);
    
    
    
    %% Seperate distinct trajectories when they merge
    
    subClusterIDs = ones(size(groupIndices));
    subTrajectoryIDs = ones(size(startIDs));
    
    % Only do this trajectories that all end up in the same place
    COHERENCE_CUTOFF = 0.7;
    bestEndCluster = mode(endClusters);
    bestStartCluster = mode(startClusters);
    if sum(bestEndCluster == endClusters) >= length(endClusters) * COHERENCE_CUTOFF && ...
            sum(bestStartCluster == startClusters) <= length(startClusters) * COHERENCE_CUTOFF
        
        %         startIndicies = find(ismember(groupIndices, startIDs));
        %         startMatrix = symmetricMatrix(startIndicies, startIndicies);
        
        startStepIDs = startIDs+1;
        startStepIDs(startStepIDs > size(asymmetricProbabilities,1)) = size(asymmetricProbabilities,1);
        startMatrix = asymmetricProbabilities(startIDs+1, startIDs+1);
        
        %         clustersLinkage = linkage(startMatrix, 'average', 'cosine');
        %
        %         rawCommunities = cluster(clustersLinkage,'maxclust',2);
        
        %         similarityMatrix = pdist2(startMatrix, startMatrix, 'cosine');
        
        rawCommunities = community_louvain(startMatrix, 0.1);
        % %         rawCommunities = GCModulMax3(startMatrix);
        %         maxRawCommunity = [];
        %         maxCost = 0;
        %         for k = 1:MAX_COMMUNITY_ATTEMPTS
        %             rawCommunities = modularity_dir(startMatrix, 0.1);
        % %             rawCommunities = community_louvain(startMatrix, 0.1);
        % %             rawCommunities = GCModulMax3(startMatrix);
        %
        %             ids = unique(rawCommunities);
        %
        %             cost = 0;
        %             innerWeight = 0;
        %             outerWeight = 0;
        %             for i = 1:length(ids)
        %                 selfIDs = rawCommunities == i;
        %
        %                 selfWeight = sum(startMatrix(selfIDs, selfIDs),2);
        %
        %                 innerWeight = innerWeight + sum(selfWeight.^2);
        %                 outerWeight = outerWeight + sum((1 - selfWeight).^2);
        %             end
        %             cost = innerWeight / outerWeight;
        %
        %             if cost > maxCost
        %                 maxCost = cost;
        %                 maxRawCommunity = rawCommunities;
        %             end
        %         end
        %
        %         rawCommunities = maxRawCommunity;
        
        bestStartClusters = [];
        coherenceLevels = [];
        clusterCounts = [];
        tempBestStartClusters = [];
        tempCoherenceLevels = [];
        tempClusterCounts = [];
        for i = 1:max(rawCommunities)
            thisClusters = startClusters(rawCommunities == i);
            tempBestStartClusters(i) = mode(thisClusters);
            tempCoherenceLevels(i) = sum(thisClusters == tempBestStartClusters(i)) / length(thisClusters);
            tempClusterCounts(i) = length(thisClusters);
        end
        
        allStarts = unique(tempBestStartClusters);
        
        for i = 1:length(allStarts)
            thisClusters = find(tempBestStartClusters == allStarts(i));
            
            bestStartClusters(i) = allStarts(i);
            clusterCounts(i) = sum(tempClusterCounts(thisClusters));
            coherenceLevels(i) = sum(tempCoherenceLevels(thisClusters) .* tempClusterCounts(thisClusters)) / clusterCounts(i);
        end
        
        [~, sortedClusters] = sort(coherenceLevels .* (clusterCounts(i) <= MIN_CLUSTER_COUNT), 'descend');
        bestStartClusters = bestStartClusters(sortedClusters);
        coherenceLevels = coherenceLevels(sortedClusters);
        clusterCounts = clusterCounts(sortedClusters);
        
        if coherenceLevels(1) >= COHERENCE_CUTOFF && clusterCounts(1) > MIN_CLUSTER_COUNT
            uniqueStarts = [];
            for i = 1:length(coherenceLevels)
                bestStart = bestStartClusters(i);
                if coherenceLevels(i) < COHERENCE_CUTOFF || clusterCounts(i) <= MIN_CLUSTER_COUNT
                    bestStarts = [];
                    for j = 1:i-1
                        bestStarts = sum(startClusters(rawCommunities == i) == bestStartClusters(j));
                    end
                    
                    [~, bestStart] = max(bestStarts);
                    bestStart = bestStartClusters(bestStart);
                    bestStartClusters(i) = bestStart;
                end
                
                uniqueStarts = unique([uniqueStarts bestStart]);
            end
            
            for i = 1:length(uniqueStarts)
                allStartClusters = find(bestStartClusters == uniqueStarts(i));
                
                allClusterValues = [];
                for j = allStartClusters
                    allClusterValues = [allClusterValues; find(rawCommunities == j)];
                end
                
                subTrajectoryIDs(allClusterValues) = i;
                
                for j = 1:length(allClusterValues)
                    thisClusterID = allClusterValues(j);
                    
                    minValue = trajectoryList{thisClusterID}(1);
                    maxValue = trajectoryList{thisClusterID}(end);
                    
                    subClusterIDs(groupIndices >= minValue & groupIndices <= maxValue) = i;
                end
            end
        end
        
    end
    
    totalSubClusters = max(subClusterIDs);
    
    if clusterID == 31
        test = 1;
    end
    
    if totalSubClusters > 1
        disp(['Splitting cluster ' num2str(clusterID) ' into ' num2str(totalSubClusters)]);
    end
    
    %%
    
    if DISPLAY_INCREMENT
        figure(5);
        clf;
        hold on;
        h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
        h.Color(4) = 0.2;
    end
    
    
    for subClusterID = 1:totalSubClusters
        %%
        subGroupIndices = find(subClusterIDs == subClusterID);
        inSubGroupIndicies = intersect(subGroupIndices, find(inGroupIndices));
        
        meanTrajectoryIndicies = subGroupIndices;
        
        subClusterPhase = localPhase(meanTrajectoryIndicies);
        subClusterInGroupPhase = localPhase(inSubGroupIndicies);
        
        subGroupStarts = find(ismember(groupIndices, startIDs(subTrajectoryIDs == subClusterID)));
        subGroupEnds = find(ismember(groupIndices, endIDs(subTrajectoryIDs == subClusterID)));
        
        %         startPhaseMean = max(-pi, min(localPhase(subGroupStarts)));
        %         endPhaseMean = min(pi, max(localPhase(subGroupEnds)));
        
        %         startPhaseMean = angle(sum(exp(1i*localPhase(subGroupStarts))));
        %         endPhaseMean = angle(sum(exp(1i*localPhase(subGroupEnds))));
        
        startPhaseMean = quantile(localPhase(subGroupStarts), 0.05);
        endPhaseMean = quantile(localPhase(subGroupEnds), 0.95);
%         startPhaseMean = mean(localPhase(subGroupStarts)) - 2*std(localPhase(subGroupStarts));
%         endPhaseMean = mean(localPhase(subGroupStarts)) + 2*std(localPhase(subGroupStarts));
        
        %         if startPhaseMean > endPhaseMean
        %             tempMean = endPhaseMean;
        %             endPhaseMean = startPhaseMean;
        %             startPhaseMean = tempMean;
        %         end
        
        SEGMENTS = ceil(abs(endPhaseMean - startPhaseMean) / USE_PHASE_VELOCITY);
        if SEGMENTS <= 0
            SEGMENTS = 1;
        end
        %         SIGMA_SCALE = 2*USE_PHASE_VELOCITY;
        sigma = SMOOTH_PHASE;
        testPhases = linspace(startPhaseMean - END_SEGMENTS*sigma, endPhaseMean + END_SEGMENTS*sigma, SEGMENTS + 2*END_SEGMENTS);
        %         testPhases = linspace(startPhaseMean, endPhaseMean, SEGMENTS);
        
        %         clusterWeights = reducedMatrix(clusterID,:);
        %         clusterWeights = clusterWeights(clusterIDs(groupIndices(meanTrajectoryIndicies)));
        
        meanValues = [];
        stdValues = [];
        weightValues = [];
        qualityValues = [];
        certaintyValues = [];
        for i = 1:SEGMENTS
            thisPhase = testPhases(i);
            offsetPhases = subClusterPhase - thisPhase;
            selfPhases = subClusterInGroupPhase - thisPhase;
            selfWeights = exp(-(selfPhases).^2/(2*sigma^2));
            
            asymmetricWeights = selfWeights'*asymmetricProbabilities(groupIndices(inSubGroupIndicies),:);
            asymmetricWeights = asymmetricWeights(groupIndices(subGroupIndices));
            %             asymmetricWeights(inSubGroupIndicies) = 1;
            
            
            weights = exp(-(offsetPhases).^2/(2*sigma^2));
            %             weights(weights < exp(-1)) = 0;
            %             weights = weights.^2;
            weights = weights .* asymmetricWeights(ismember(subGroupIndices, meanTrajectoryIndicies))';
            weights = weights .* stateValidities(groupIndices(meanTrajectoryIndicies));
            %             if sum(weights) > 0
            weights = weights / sum(weights);
            %             end
            
            validIndicies = 1:length(meanTrajectoryIndicies);
            
            %             weightedSTD = sqrt(var(finalDynamicsStream(groupIndices(subGroupIndices),:), weights));
            %             weightedStream = finalDynamicsStream(groupIndices(subGroupIndices),:) ./ weightedSTD;
            %             weightedMean = sum(weights.*weightedStream);
            %             weightedDistances = sqrt(sum((weightedStream - weightedMean).^2, 2));
            %             weightedDistances = weightedDistances(weights > 1 / length(weights));
            %
            %             smoothSigma = length(weightedDistances)/10;
            %             x = -ceil(3*smoothSigma):ceil(3*smoothSigma);
            %             kernel = -(smoothSigma^2 - x.^2).*exp(-x.^2/(2*smoothSigma^2))/(smoothSigma^5*sqrt(2*pi));
            % %             kernel = -x.*exp(-x.^2/(2*smoothSigma^2))/(smoothSigma^3*sqrt(2*pi));
            % %             kernel = exp(-x.^2/(2*smoothSigma^2));
            % %             kernel = kernel / sum(kernel);
            %
            % %             smoothedDistances = filterData(sort(weightedDistances'), smoothSigma, [], 1, 0);
            %             [sortedWeights, sortIDs] = sort(weightedDistances');
            %             smoothedCurvature = filterData(sortedWeights, 0, kernel, 1, 0);
            %
            %             fitMatrix = [ones(length(weightedDistances), 1)'; 1:length(weightedDistances)]';
            %             b = fitMatrix \ sortedWeights';
            %             bestFitLine = b(1)+b(2)*(1:length(weightedDistances));
            %
            %             residues = sort(weightedDistances)' - bestFitLine;
            %             residueSTD = std(residues);
            %
            %             firstOutlier = find(abs(residues) <= residueSTD*2, 1, 'last');
            %             if isempty(firstOutlier)
            %                 firstOutlier = 0;
            %             end
            %
            %             %If there are a LOT of outliers, it's probably not actually an
            %             %outlier...
            %             if firstOutlier > length(weightedDistances) * 0.8
            % %                 disp([num2str(i) ') Removing ' num2str(length(firstOutlier:length(sortIDs))) ' outliers']);
            %                 validIndicies(sortIDs(firstOutlier:length(sortIDs))) = [];
            %             end
            %
            % %             figure(1);
            % %             clf;
            % %             hold on
            % %             plot(sort(weightedDistances'));
            % %             plot(bestFitLine);
            
            if sum(weights) == 0
                test = 2;
            end
            
            meanValues(i,:) = sum(weights(validIndicies).*finalDynamicsStream(groupIndices(meanTrajectoryIndicies(validIndicies)),:));
            stdValues(i,:) = sqrt(var(finalDynamicsStream(groupIndices(meanTrajectoryIndicies(validIndicies)),:), weights(validIndicies))) / sqrt(sum(weights(validIndicies) > 1 / length(validIndicies)));
            weightValues(i,:) = weights;
            inClusterWeight = sum(weights(ismember(subGroupIndices, inSubGroupIndicies)));
            outClusterWeight = sum(weights(ismember(subGroupIndices, setdiff(subGroupIndices, inSubGroupIndicies))));
            qualityValues(i) = inClusterWeight / (inClusterWeight + outClusterWeight);
            certaintyValues(i) = mean(weights(ismember(subGroupIndices, inSubGroupIndicies)));
            
            %             scatter3(finalDynamicsStream(groupIndices(meanTrajectoryIndicies),:)*pcaBasis(:,1), finalDynamicsStream(groupIndices(meanTrajectoryIndicies),:)*pcaBasis(:,2), finalDynamicsStream(groupIndices(meanTrajectoryIndicies),:)*pcaBasis(:,3), 32, weights);
            % %             scatter3(finalDynamicsStream(groupIndices(meanTrajectoryIndicies),:)*pcaBasis(:,1), finalDynamicsStream(groupIndices(meanTrajectoryIndicies),:)*pcaBasis(:,2), finalDynamicsStream(groupIndices(meanTrajectoryIndicies),:)*pcaBasis(:,3), 32, asymmetricWeights);
            %             scatter3(meanValues(end,:)*pcaBasis(:,1), meanValues(end,:)*pcaBasis(:,2), meanValues(end,:)*pcaBasis(:,3), 32, 'kx');
            
            %             if i == 1 && subClusterID == 2
            %                 test = 1;
            %             end
        end
        
        finalPhases = nan(size(finalDynamicsStream,1),1);
        finalPhases(groupIndices(subGroupIndices)) = subClusterPhase;
        
        trajectoryClusterMap(end+1) = clusterID;
        trajectoryMeans{length(trajectoryClusterMap)} = meanValues;
        trajectorySTDs{length(trajectoryClusterMap)} = stdValues;
        trajectoryQuality{end+1} = qualityValues;
        trajectoryCertainty{end+1} = certaintyValues;
        trajectoryCounts(end+1) = length(trajectoryList);
        trajectoryIndices{end+1} = groupIndices(inSubGroupIndicies);
        trajectoryPhases(:, length(trajectoryClusterMap)) = finalPhases;
        dataTrajectoryMap(groupIndices(inSubGroupIndicies)) = length(trajectoryClusterMap);
        
        allWeights = zeros(size(finalDynamicsStream,1), size(weightValues,1));
        allWeights(groupIndices(subGroupIndices), :) = weightValues';
        trajectoryWeights{length(trajectoryClusterMap)} = allWeights;
        
        if DISPLAY_INCREMENT
            %         scatter3(finalDynamicsStream(groupIndices(subGroupIndices),:)*pcaBasis(:,1), finalDynamicsStream(groupIndices(subGroupIndices),:)*pcaBasis(:,2), finalDynamicsStream(groupIndices(subGroupIndices),:)*pcaBasis(:,3), 32, localPhase(subGroupIndices));
            scatter3(finalDynamicsStream(groupIndices(meanTrajectoryIndicies),:)*pcaBasis(:,1), finalDynamicsStream(groupIndices(meanTrajectoryIndicies),:)*pcaBasis(:,2), finalDynamicsStream(groupIndices(meanTrajectoryIndicies),:)*pcaBasis(:,3), 32, localPhase(meanTrajectoryIndicies));
            plot3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 'LineWidth', 3);
            scatter3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 32, testPhases);
        end
    end
    
    if length(trajectoryMeans) == 53
        tersdt = 1;
    end
    
    if DISPLAY_INCREMENT
        title(num2str(clusterID));
        
        %         caxis([-pi pi]);
        drawnow
        pause(1);
    else
        waitHandle.iterate(1);
    end
end
if ~DISPLAY_INCREMENT
    close(waitHandle);
end

%%

if ~DISPLAY_INCREMENT
    figure(4);
    clf;
    hold on;
    h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
    h.Color(4) = 0.2;
    for i = 1:length(trajectoryMeans)
        meanValues = trajectoryMeans{i};
        %         plot3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 'LineWidth', 3);
        %         scatter3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 32, trajectoryQuality(:,i).^2);
        scatter3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 32, trajectoryCertainty{i});
        %         scatter3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 32, ones(size(trajectoryCertainty{i}))*maxSelfTransition(i));
    end
    colormap(jet(256))
end

%%

% Get detailed balance decomp
[steadyState, ~] = eigs(asymmetricProbabilities', 1);
fluxMatrix = diag(steadyState) * asymmetricProbabilities;
symmetricMatrix = (fluxMatrix+fluxMatrix.')/2;
antisymmetricMatrix = (fluxMatrix-fluxMatrix.')/2;
symmetricMatrix = symmetricMatrix ./ steadyState;
antisymmetricMatrix = antisymmetricMatrix ./ steadyState;

% softAsymmetricProbabilities = asymmetricProbabilities^3;

stateTransitionMap = [];
stateWeights = [];
stateSimilarities = [];
stateTrajectoryMap = [];
trajectoryStateMap = {};
stateID = 1;
waitHandle = parfor_progressbar(length(trajectoryWeights), 'Building trajectory transition map');
for i = 1:length(trajectoryWeights)
    thisWeights = trajectoryWeights{i};
    
    for j = 1:size(thisWeights,2)
        stateTrajectoryMap(stateID,:) = [i, j];
        trajectoryStateMap{i}(j) = stateID;
        
        selfWeights = thisWeights(:,j);
        asymmetricWeights = selfWeights'*asymmetricProbabilities;
        asymmetricWeights = asymmetricWeights / sum(asymmetricWeights);
        
        if sum(isnan(asymmetricWeights(:))) > 0
            test = 1;
        end
        if sum(isnan(selfWeights(:))) > 0
            test = 1;
        end
        
        stateWeights(stateID, :) = selfWeights;
        stateTransitionMap(stateID, :) = asymmetricWeights;
        
        stateID = stateID + 1;
    end
    
    waitHandle.iterate(1);
end
close(waitHandle);


stateSimilarities = stateTransitionMap * stateWeights';
stateSimilarities = stateSimilarities ./ sum(stateSimilarities,2);

% asymmetricStateSimilarities = stateSimilarities;
% asymmetricStateSimilarities(1,:) = [];
% asymmetricStateSimilarities(:,end) = [];

% empricalStateTransitionMap = [];
% for i = 1:size(stateTrajectoryMap,1)
%
% end

clf;
imagesc(stateSimilarities);

%% Find reduced phases

% Get detailed balance decomp
[steadyState, ~] = eigs(stateSimilarities', 1);
fluxMatrix = diag(steadyState) * stateSimilarities;
symmetricMatrix = (fluxMatrix+fluxMatrix.')/2;
antisymmetricMatrix = (fluxMatrix-fluxMatrix.')/2;
symmetricMatrix = symmetricMatrix ./ steadyState;
antisymmetricMatrix = antisymmetricMatrix ./ steadyState;

% Rescale diffusion map
diffusedProbabilities = max(0, symmetricMatrix + sign(antisymmetricMatrix).*abs(antisymmetricMatrix).^(1/4));
diffusedProbabilities = diffusedProbabilities ./ sum(diffusedProbabilities,2);

[eigenVectors, eigenValues] = eigs(diffusedProbabilities, 10);
eigenValues = diag(eigenValues);

[~, bestComplexEigenVector] = max((abs(imag(eigenValues)) > 0) .* abs(eigenValues));


reducedGobalPhase = angle(eigenVectors(:,bestComplexEigenVector));


figure(5);
clf;
hold on;
h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
h.Color(4) = 0.2;
for i = 1:length(trajectoryMeans)
    meanValues = trajectoryMeans{i};
    stateIDs = trajectoryStateMap{i};
    
    scatter3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 32, reducedGobalPhase(stateIDs));
end
colormap(hsv);

%%

% startConnections = stateSimilarities(1:size(trajectoryWeights,2):end,:);
% endConnections = stateSimilarities(30:size(trajectoryWeights,2):end,:);
% for i = 1:size(startConnections,1)
%     for j = 1:size(startConnections,1)
%         thisIDs = trajectoryStateMap(j,:);
%
%         startConnections(i,j) = sum(startConnections(i,thisIDs));
%         endConnections(i,j) = sum(endConnections(i,thisIDs));
%     end
% end
% startConnections(:,size(startConnections,1)+1:end) = [];
% endConnections(:,size(endConnections,1)+1:end) = [];

%%

END_TIME = 2;

allEndTrajectories = {};
allStartTrajectories = {};
allEndTrajectoryWeights = {};
allStartTrajectoryWeights = {};
maxSelfTransition = [];
startConnections = [];
endConnections = [];
startConnectionsWeighted = [];
endConnectionsWeighted = [];
for i = 1:length(trajectoryMeans)
    maxSelfTransition(i) = max(trajectoryQuality{i});
    
    % end weights
    
    endWeights = trajectoryWeights{i}(:,end-END_SEGMENTS);
    endWeights = [zeros(END_TIME,1); endWeights];
    endWeights(end-END_TIME+1:end) = [];
    endWeights(trajectoryIndices{i}) = 0;
    
    endStateConnections = endWeights' * stateWeights';
    for j = 1:length(trajectoryMeans)
        thisIDs = trajectoryStateMap{j};
        
        endStateConnections(j) = sum(endStateConnections(thisIDs)) .* trajectoryCounts(j)^2;
    end
    endStateConnections(length(trajectoryMeans)+1:end) = [];
    endStateConnections = endStateConnections /sum(endStateConnections);
    allEndTrajectories{i} = find(endStateConnections > 0.001);
    allEndTrajectoryWeights{i} = endStateConnections(allEndTrajectories{i});
    endConnectionsWeighted(i,:) = endStateConnections;
    
    thisIndicies = trajectoryIndices{i} + 1;
    thisIndicies = setdiff(thisIndicies, trajectoryIndices{i});
    thisIndicies(thisIndicies > size(dataTrajectoryMap,2)) = [];
    connections = tabulate(dataTrajectoryMap(thisIndicies));
    endStateConnections = zeros(1, length(trajectoryMeans));
    endStateConnections(1:size(connections,1)) = connections(:,3) / 100;
    endConnections(i,:) = endStateConnections;
    
    
    %     checkIDs = find(endWeights > 0.001);
    %
    %     endTrajectories = [];
    %     endTrajectoryWeights = [];
    %     for j = 1:length(checkIDs)
    %         thisID = checkIDs(j);
    %         bestTrajectoryID = find(trajectoryClusterMap == clusterIDs(thisID));
    %
    %         if isempty(bestTrajectoryID)
    %             continue;
    %         end
    %
    %         trajectoryID = find(endTrajectories == bestTrajectoryID);
    %         if isempty(trajectoryID)
    %             trajectoryID = length(endTrajectories) + 1;
    %
    %             endTrajectories(trajectoryID) = bestTrajectoryID;
    %             endTrajectoryWeights(trajectoryID) = 0;
    %         end
    %
    %         endTrajectoryWeights(trajectoryID) = endTrajectoryWeights(trajectoryID) + endWeights(thisID);
    %     end
    %     endTrajectoryWeights = endTrajectoryWeights / sum(endTrajectoryWeights);
    %
    %     allEndTrajectories{i} = endTrajectories;
    %     allEndTrajectoryWeights{i} = endTrajectoryWeights;
    
    % start weights
    
    startWeights = trajectoryWeights{i}(:,1+END_SEGMENTS);
    startWeights = [startWeights; zeros(END_TIME,1)];
    startWeights(1:END_TIME) = [];
    startWeights(trajectoryIndices{i}) = 0;
    
    startStateConnections = startWeights' * stateWeights';
    for j = 1:length(trajectoryMeans)
        thisIDs = trajectoryStateMap{j};
        
        startStateConnections(j) = sum(startStateConnections(thisIDs)) .* trajectoryCounts(j)^2;
    end
    startStateConnections(length(trajectoryMeans)+1:end) = [];
    startStateConnections = startStateConnections /sum(startStateConnections);
    allStartTrajectories{i} = find(startStateConnections > 0.001);
    allStartTrajectoryWeights{i} = startStateConnections(allStartTrajectories{i});
    startConnectionsWeighted(i,:) = startStateConnections;
    
    thisIndicies = trajectoryIndices{i} - 1;
    thisIndicies = setdiff(thisIndicies, trajectoryIndices{i});
    thisIndicies(thisIndicies < 1) = [];
    connections = tabulate(dataTrajectoryMap(thisIndicies));
    startStateConnections = zeros(1, length(trajectoryMeans));
    startStateConnections(1:size(connections,1)) = connections(:,3) / 100;
    startConnections(i,:) = startStateConnections;
    
    %     checkIDs = find(startWeights > 0.001);
    %
    %     startTrajectories = [];
    %     startTrajectoryWeights = [];
    %     for j = 1:length(checkIDs)
    %         thisID = checkIDs(j);
    %         bestTrajectoryID = find(trajectoryClusterMap == clusterIDs(thisID));
    %
    %         if isempty(bestTrajectoryID)
    %             continue;
    %         end
    %
    %         trajectoryID = find(startTrajectories == bestTrajectoryID);
    %         if isempty(trajectoryID)
    %             trajectoryID = length(startTrajectories) + 1;
    %
    %             startTrajectories(trajectoryID) = bestTrajectoryID;
    %             startTrajectoryWeights(trajectoryID) = 0;
    %         end
    %
    %         startTrajectoryWeights(trajectoryID) = startTrajectoryWeights(trajectoryID) + startWeights(thisID);
    %     end
    %     startTrajectoryWeights = startTrajectoryWeights / sum(startTrajectoryWeights);
    %
    %     allStartTrajectories{i} = startTrajectories;
    %     allStartTrajectoryWeights{i} = startTrajectoryWeights;
    %
    %     thisStartConnections = zeros(1, length(trajectoryMeans));
    %     thisStartConnections(allStartTrajectories{i}) = allStartTrajectoryWeights{i};
    %     startConnections(i,:) = thisStartConnections;
    %
    %     thisEndConnections = zeros(1, length(trajectoryMeans));
    %     thisEndConnections(allEndTrajectories{i}) = allEndTrajectoryWeights{i};
    %     endConnections(i,:) = thisEndConnections;
    
    %     figure(4);
    %     clf;
    %     hold on;
    %     colormap(jet(256));
    %     h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
    %     h.Color(4) = 0.2;
    %     plotPoints = find(endWeights > 0.001);
    %     scatter3(finalDynamicsStream(plotPoints,:)*pcaBasis(:,1), finalDynamicsStream(plotPoints,:)*pcaBasis(:,2), finalDynamicsStream(plotPoints,:)*pcaBasis(:,3), 32, endWeights(plotPoints));
    
end

transposedEnds = endConnectionsWeighted';
transposedEnds = transposedEnds ./ sum(transposedEnds,2);
mergedStartConnections = transposedEnds .* startConnectionsWeighted;
mergedStartConnections = mergedStartConnections .* ~eye(size(mergedStartConnections));
transposedStarts = startConnectionsWeighted';
transposedStarts = transposedStarts ./ sum(transposedStarts,2);
mergedEndConnections = endConnectionsWeighted .* transposedStarts;
mergedEndConnections = mergedEndConnections .* ~eye(size(mergedEndConnections));

% figure(2);
% clf;
% colormap(jet(256));
% subplot(2,2,1);
% imagesc(endConnections);
% title("End connections");
% subplot(2,2,2);
% imagesc(startConnections);
% title("Start connections");
% % subplot(2,2,3);
% % imagesc(endConnectionsWeighted);
% % title("Merge start connections");
% % subplot(2,2,4);
% % imagesc(startConnectionsWeighted);
% % title("Merge end connections");
% subplot(2,2,3);
% imagesc(mergedStartConnections);
% title("Merge start connections");
% subplot(2,2,4);
% imagesc(mergedEndConnections);
% title("Merge end connections");
% 
% figure(4);
% clf;
% hold on;
% h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
% h.Color(4) = 0.2;
% for i = 1:length(trajectoryMeans)
%     meanValues = trajectoryMeans{i};
%     %         plot3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 'LineWidth', 3);
%     %         scatter3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 32, trajectoryQuality(:,i).^2);
%     %         scatter3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 32, trajectoryCertainty(:,i));
%     %     scatter3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 32, ones(size(trajectoryCertainty(:,i)))*length(allStartTrajectories{i}));
%     scatter3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 32, ones(size(trajectoryCertainty{i}))*length(allEndTrajectories{i}));
% end




%%


% t = minspantree(g, 'Type', 'forest');
% nonTreeEdges = setdiff(g.Edges.EndNodes, t.Edges.EndNodes, 'rows');
% cycles = cell(size(nonTreeEdges, 1), 1);
% cycleWeights = [];
% for ii=1:length(cycles)
%     src = nonTreeEdges(ii, 1); tgt = nonTreeEdges(ii, 2);
%     [path, pathWeight, edgepath] = shortestpath(t, src, tgt);
%     cycles{ii} = [tgt path];
%
%     cycleWeights(ii) = pathWeight / length(path);
% end
%
% [~, cycleOrder] = sort(cycleWeights);

%%

% PLOT_CYCLE = 2;
%
% figure(5);
% clf;
% hold on;
% h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
% h.Color(4) = 0.2;
% cyclePath = cycles{cycleOrder(PLOT_CYCLE)};
% for i = 1:length(cyclePath) - 1
%     plot3(trajectoryMeans(:,:,cyclePath(i))*pcaBasis(:,1), trajectoryMeans(:,:,cyclePath(i))*pcaBasis(:,2), trajectoryMeans(:,:,cyclePath(i))*pcaBasis(:,3), 'LineWidth', 3);
%
%     points(1,:) = trajectoryMeans(end,:,cyclePath(i));
%     points(2,:) = trajectoryMeans(1,:,cyclePath(i+1));
%
%     plot3(points*pcaBasis(:,1), points*pcaBasis(:,2), points*pcaBasis(:,3), 'k', 'LineWidth', 1);
% end

%%

CONNECT_CUTOFF = 0.01;

normWeights = stateWeights ./ sqrt(sum(stateWeights.^2,2));
weightOverlap = normWeights * normWeights';

stateOverlap = [];
crossStateOverlap = [];
for i = 1:length(trajectoryMeans)
    thisIDs = trajectoryStateMap{i};
    
    for j = 1:length(trajectoryMeans)
        thoseIDs = trajectoryStateMap{j};
        stateOverlap(i,j) = sum(sum(weightOverlap(thisIDs,thoseIDs)));
    end
    stateOverlap(i,:) = stateOverlap(i,:) / sum(stateOverlap(i,:));
    
    thisOverlap = stateOverlap(i,:);% / stateOverlap(i,i);
    [~, sortOrder] = sort(thisOverlap, 'descend');
    %     thisOverlap(sortOrder(1)) = 0;
    %     nonSelf = sum(thisOverlap);
    thisOverlap(sortOrder(1:3)) = 0;
    crossStateOverlap(i) = sum(thisOverlap);
end


figure(5);
clf;
hold on;
h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
h.Color(4) = 0.2;
colorValues = log(crossStateOverlap);
% colorValues(colorValues < -10) = -10;
caxis([min(colorValues) max(colorValues)]);
colorValues = colorValues - min(colorValues);
colorValues = colorValues / max(colorValues);
colorValues = round(colorValues*255 + 1);
colors = jet(256);
colormap(jet(256));
colorbar
for i = 1:length(trajectoryMeans)
    %     if sum(mergedTrajectories(i,:)) + sum(mergedTrajectories(:,i)) <= 0
    %         continue;
    %     end
    
    
    maxValue = max(mergedStartConnections(i,:));
    startTrajectories = find(mergedStartConnections(i,:) / maxValue > CONNECT_CUTOFF);% & mergedTrajectories(i,:));
    
    maxValue = max(mergedEndConnections(i,:));
    endTrajectories = find(mergedEndConnections(i,:) / maxValue > CONNECT_CUTOFF);% & mergedTrajectories(i,:));
    
    for j = 1:length(startTrajectories)
        points(1,:) = trajectoryMeans{i}(1,:);
        points(2,:) = trajectoryMeans{startTrajectories(j)}(end,:);
        
        plot3(points*pcaBasis(:,1), points*pcaBasis(:,2), points*pcaBasis(:,3), 'k', 'LineWidth', 1);
    end
    
    for j = 1:length(endTrajectories)
        points(1,:) = trajectoryMeans{i}(end,:);
        points(2,:) = trajectoryMeans{endTrajectories(j)}(1);
        
        plot3(points*pcaBasis(:,1), points*pcaBasis(:,2), points*pcaBasis(:,3), 'k', 'LineWidth', 1);
    end
    
    plot3(trajectoryMeans{i}*pcaBasis(:,1), trajectoryMeans{i}*pcaBasis(:,2), trajectoryMeans{i}*pcaBasis(:,3), 'Color', colors(colorValues(i),:), 'LineWidth', 3);
    
end

% %% Find best phase alignemnts
%
% bestPhaseOffsets = [];
% for i = 1:size(trajectoryPhases, 2)
%     thisPhases = trajectoryPhases(:,i);
%
%     for j = 1:size(trajectoryPhases, 2)
%         comparePhases = trajectoryPhases(:,j);
%
%         phaseFit = @(x) nansum(angleDiff(comparePhases, (thisPhases + x)).^2);
%
%         options = optimoptions('fmincon','Display','off');
%         bestPhaseOffset = fmincon(phaseFit, 0, [], [], [], [], -pi, pi, [], options);
%
%         bestPhaseOffsets(i,j) = bestPhaseOffset;
%     end
% end
%
% TEST_TRAJECTORY = 10;
% figure(5);
% clf;
% hold on;
% h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
% h.Color(4) = 0.2;
% for i = 1:length(trajectoryIndices)
%     thisPoints = trajectoryIndices{i};
%
%     if stateOverlap(TEST_TRAJECTORY,i) > 0.01
%         scatter3(finalDynamicsStream(thisPoints,:)*pcaBasis(:,1), finalDynamicsStream(thisPoints,:)*pcaBasis(:,2), finalDynamicsStream(thisPoints,:)*pcaBasis(:,3), 32, wrapToPi(trajectoryPhases(thisPoints,i) - bestPhaseOffsets(TEST_TRAJECTORY,i)));
%     end
% end
% colormap(hsv(256))


%% Find loops greedy

% startTrajectoryCount = startConnections ./ max(startConnections,[],2) .* trajectoryCounts';
% endTrajectoryCount = endConnections ./ max(endConnections,[],2).* trajectoryCounts';

% connectionTrajectoryCount = ;
% bestConnections = bestConnections ./ max(bestConnections,[],2);

% graph = grGraphFromMatrix(connectionTrajectoryCount);

% minTree = grCycleBasis(graph(:,1:2));

MIN_VALUE = 1/10000;

graphTree = (endConnectionsWeighted);
graphTree(graphTree < MIN_VALUE) = MIN_VALUE;
graphTree(find(eye(size(graphTree)))) = MIN_VALUE;
graphTree = 1 ./ graphTree;

% g = digraph(graphTree);

if DISPLAY
    figure(5);
    clf;
    hold on;
    h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
    h.Color(4) = 0.2;
end

allLoops  = {};
loopCorrelationMatrix = zeros(size(startConnectionsWeighted, 1));
loopValues = [];
loopCounts = zeros(1, size(startConnectionsWeighted, 1));
for i = 1:size(startConnectionsWeighted,1)
    loopStart = i;
    possibleEnds = find(startConnectionsWeighted(loopStart,:));
    bestLength = 99999999;
    bestPath = [];
    for j = 1:length(possibleEnds)
        if possibleEnds(j) == loopStart
            continue;
        end
        
        [pathLength, path, pred] = graphshortestpath(sparse(graphTree), loopStart, possibleEnds(j));
        
        totalLength = pathLength + graphTree(possibleEnds(j), loopStart);
        
        if totalLength < bestLength
            bestLength = totalLength;
            bestPath = path;
        end
    end
    
    allLoops{loopStart} = bestPath;
    loopCorrelationMatrix(loopStart,bestPath) = 1;
    loopValues(loopStart) = bestLength;
    loopCounts(bestPath) = loopCounts(bestPath) + 1;
    
    if DISPLAY
        drawPath = [bestPath bestPath(1)];
        for j = 1:length(drawPath)-1
            thisID = drawPath(j);
            
            points(1,:) = trajectoryMeans{thisID}(end,:);
            points(2,:) = trajectoryMeans{drawPath(j+1)}(1,:);
            plot3(points*pcaBasis(:,1), points*pcaBasis(:,2), points*pcaBasis(:,3), 'k', 'LineWidth', 2);
            
            plot3(trajectoryMeans{thisID}*pcaBasis(:,1), trajectoryMeans{thisID}*pcaBasis(:,2), trajectoryMeans{thisID}*pcaBasis(:,3), 'LineWidth', 3);
        end
    end
    
end

loopOverlapMatrix = stateOverlap ./ max(stateOverlap, [], 2);

finalLoopOverlaps = [];
for i = 1:size(loopOverlapMatrix,1)
    thisStates = find(loopCorrelationMatrix(i,:));
    for j = 1:size(loopOverlapMatrix,1)
        otherStates = find(loopCorrelationMatrix(j,:));
        
        overlapStates = intersect(thisStates, otherStates);
        
        if length(overlapStates) > 0
            finalLoopOverlaps(i,j) = length(overlapStates);
        else
            finalLoopOverlaps(i,j) = 0;
        end
    end
end
finalLoopOverlaps = 1./finalLoopOverlaps;
finalLoopOverlaps = finalLoopOverlaps .* ~eye(size(finalLoopOverlaps));
finalLoopOverlaps = finalLoopOverlaps ./ diag(finalLoopOverlaps);
finalLoopOverlaps = 1 - finalLoopOverlaps;
finalLoopOverlaps(finalLoopOverlaps<0) = 0;
finalLoopOverlaps(isnan(finalLoopOverlaps)) = 0;

% loopClusteringMatrix = pdist2(loopCorrelationMatrix,loopCorrelationMatrix, 'cityblock');

loopLinkage = linkage(squareform(finalLoopOverlaps), 'average' );
loopCommunities = cluster(loopLinkage,'maxclust',8);
% loopCommunities = link_communities(loopCorrelationMatrix, 'single' );
% loopCommunities = community_louvain(loopCorrelationMatrix, 0.9);

% [~,sortOrder] = sort(loopCommunities);
% clf
% imagesc(loopCorrelationMatrix(sortOrder,sortOrder));


% Calculate loop weights
totalTrajectoryCount = 0;

loopWeights = [];
for i = 1:max(loopCommunities)
    thisLoopIDs = find(loopCommunities == i);
    
    totalTrajectories = [];
    
    notLoopIDs = 1:size(startConnectionsWeighted, 1);
    inLoopIDs = zeros(1, size(startConnectionsWeighted, 1));
    thisWeights = zeros(1, size(startConnectionsWeighted, 1));
    for j = 1:length(thisLoopIDs)
        loopSubtrajectoryID = thisLoopIDs(j);
        
        subLoops = allLoops{loopSubtrajectoryID};
        
        notLoopIDs(subLoops) = 0;
        inLoopIDs(subLoops) = 1;
        
        thisWeights(subLoops) = thisWeights(subLoops) + 1;
        
        totalTrajectories(j) = sum(trajectoryCounts(subLoops));
        
        totalTrajectoryCount = totalTrajectoryCount + totalTrajectories(j);
    end
    thisWeights = thisWeights / sum(thisWeights);
    loopWeights(i,:) = thisWeights;
    
    inLoopIDs = find(inLoopIDs);
end

PLOT_LOOPS = 1:max(loopCommunities);

% Get loop point weights
loopPointWeights = zeros(length(PLOT_LOOPS), size(asymmetricProbabilities,1));
for i = 1:max(loopCommunities)
    subLoops = find(loopWeights(i,:) > 0);
    
    allStates = [];
    for j = 1:length(subLoops)
        thisTrajectory = subLoops(j);
        thisStates = trajectoryStateMap{thisTrajectory};
        
        allStates = [allStates, thisStates];
    end
    
    loopPointWeight = zeros(1, size(asymmetricProbabilities,1));
    for j = 1:length(allStates)
        trajectoryID = stateTrajectoryMap(allStates(j), 1);
        
        trajectoryWeight = loopWeights(i, trajectoryID);
        
        maxWeights = stateWeights(allStates(j),:);
        maxWeights = maxWeights / max(maxWeights);
        
        loopPointWeight = loopPointWeight + maxWeights * trajectoryWeight;
    end
    
    thisLoopPoints = find(loopPointWeight > 0.001);
    loopPointWeights(i,thisLoopPoints) = loopPointWeight(thisLoopPoints);
    
    
    figure(i);
    clf;
    hold on;
    h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
    h.Color(4) = 0.2;
    %     scatter3(finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,1), finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,2), finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,3), 32, currentDistribution);
        scatter3(finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,1), finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,2), finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,3), 32, loopPointWeight(thisLoopPoints));
    colormap(parula(256));
    subLoops = find(loopWeights(i,:) > 0);
    
    for j = subLoops
        thisStates = trajectoryStateMap{j};
    end
    
end

% %%
% 
% PLOT_LOOPS = 1:max(loopCommunities);
% % PLOT_LOOPS = 4;




%% Calculate loop weights
totalTrajectoryCount = 0;

loopCommunities = 1:size(stateOverlap,1);

loopWeights = [];
for i = 1:max(loopCommunities)
    thisLoopIDs = find(loopCommunities == i);
    
    totalTrajectories = [];
    
    notLoopIDs = 1:size(startConnectionsWeighted, 1);
    inLoopIDs = zeros(1, size(startConnectionsWeighted, 1));
    thisWeights = zeros(1, size(startConnectionsWeighted, 1));
    for j = 1:length(thisLoopIDs)
        loopSubtrajectoryID = thisLoopIDs(j);
        
        subLoops = allLoops{loopSubtrajectoryID};
        
        notLoopIDs(subLoops) = 0;
        inLoopIDs(subLoops) = 1;
        
        thisWeights(subLoops) = thisWeights(subLoops) + 1;
        
        totalTrajectories(j) = sum(trajectoryCounts(subLoops));
        
        totalTrajectoryCount = totalTrajectoryCount + totalTrajectories(j);
    end
    thisWeights = thisWeights / sum(thisWeights);
    loopWeights(i,:) = thisWeights;
    
    inLoopIDs = find(inLoopIDs);
end

%% Get all loop qualities
USE_ALL_LOOPS = 1;
MIN_OVERLAP = 0.05;

loopCommunities = 1:size(stateOverlap,1);

totalTrajectoryCount = 0;
loopQualities = [];
loopCertainties = [];
loopOverlap = [];
for i = 1:max(loopCommunities)
    thisLoopIDs = find(loopCommunities == i);
    
    totalTrajectories = [];
    allPoints = [];
    allStates = [];
    
    notLoopIDs = 1:size(startConnectionsWeighted, 1);
    inLoopIDs = zeros(1, size(startConnectionsWeighted, 1));
    for j = 1:length(thisLoopIDs)
        loopSubtrajectoryID = thisLoopIDs(j);
        
        subLoops = allLoops{loopSubtrajectoryID};
        
        notLoopIDs(subLoops) = 0;
        inLoopIDs(subLoops) = 1;
        
        totalTrajectories(j) = sum(trajectoryCounts(subLoops));
        
        totalTrajectoryCount = totalTrajectoryCount + totalTrajectories(j);
        
        if ~USE_ALL_LOOPS
            allPoints = unique([allPoints trajectoryIndices{loopSubtrajectoryID}']);
            allStates = unique([allStates trajectoryStateMap{loopSubtrajectoryID}]);
        else
            for k = subLoops
                allPoints = unique([allPoints trajectoryIndices{k}']);
                allStates = unique([allStates trajectoryStateMap{k}]);
            end
        end
    end
    
    if USE_ALL_LOOPS
        inLoopIDs = find(inLoopIDs);% Use ALL loop trajectories in cluster
    else
        inLoopIDs = thisLoopIDs; % Use only trajectories in cluster
    end
    
    %%
    
    notLoopIDs(notLoopIDs == 0) = [];
    thisOverlaps = zeros(1, size(startConnectionsWeighted, 1));
    for j = 1:size(startConnectionsWeighted, 1)
        thisOverlaps(j) = sum(stateOverlap(j, notLoopIDs));
    end
    
    loopCertainties(i) = mean(totalTrajectories);
    loopOverlap(i) = max(MIN_OVERLAP, sum(loopWeights(i,:) .* (thisOverlaps)));
    loopQualities(i) = loopCertainties(i) / loopOverlap(i);
end

if length(loopQualities) > 1
    [loopQualityClusters, centroids] = kmeans(log(loopQualities'), 2);
    
    bestThreshold = (exp(mean(centroids)));
    [~, bestCentroid] = max(centroids);
else
    loopQualityClusters = 1;
    bestCentroid = 1;
end

% goodLoops = find(loopQualityClusters == bestCentroid);
% finalCommunities = loopCommunities;

%% Find best loop count
MIN_VALUE = 1/10000;

finalLoopOverlaps = [];
for i = 1:size(stateOverlap,1)
    thisStates = allLoops{i};
    for j = 1:size(stateOverlap,1)
        otherStates = allLoops{j};
        
        overlapStates = intersect(thisStates, otherStates);
        
        if length(overlapStates) > 0
            finalLoopOverlaps(i,j) = min(max(stateOverlap(thisStates, otherStates), [], 2));
        else
            finalLoopOverlaps(i,j) = 0;
        end
    end
end
finalLoopOverlaps(finalLoopOverlaps < MIN_VALUE) = MIN_VALUE;

% loopOverlapCounts = pdist2(loopCorrelationMatrix,loopCorrelationMatrix, 'jaccard') .* loopQualities';
% loopOverlapCounts = pdist2(loopCorrelationMatrix,loopCorrelationMatrix, 'jaccard') .* loopQualities';
loopOverlapCounts = 1 ./ finalLoopOverlaps;
loopOverlapCounts = loopOverlapCounts .* ~eye(size(loopOverlapCounts));
loopOverlapCounts = loopOverlapCounts .* loopQualities';
loopOverlapCounts = min(loopOverlapCounts, loopOverlapCounts');

loopLinkage = linkage(squareform(loopOverlapCounts), 'complete' );

bestClusterCount = 0;
bestLoopQuality = 0;
loopTotalQualities = [];

for clusterCounts = 1:10
    loopCommunities = cluster(loopLinkage,'maxclust',clusterCounts);

    %% Calculate loop weights
    totalTrajectoryCount = 0;

    loopWeights = [];
    for i = 1:max(loopCommunities)
        thisLoopIDs = find(loopCommunities == i);

        totalTrajectories = [];

        notLoopIDs = 1:size(startConnectionsWeighted, 1);
        inLoopIDs = zeros(1, size(startConnectionsWeighted, 1));
        thisWeights = zeros(1, size(startConnectionsWeighted, 1));
        for j = 1:length(thisLoopIDs)
            loopSubtrajectoryID = thisLoopIDs(j);

            subLoops = allLoops{loopSubtrajectoryID};

            notLoopIDs(subLoops) = 0;
            inLoopIDs(subLoops) = 1;

            thisWeights(subLoops) = thisWeights(subLoops) + 1;

            totalTrajectories(j) = sum(trajectoryCounts(subLoops));

            totalTrajectoryCount = totalTrajectoryCount + totalTrajectories(j);
        end
        thisWeights = thisWeights / sum(thisWeights);
        loopWeights(i,:) = thisWeights;

        inLoopIDs = find(inLoopIDs);
    end

    %% Get loop quality
    USE_ALL_LOOPS = 1;
    MIN_OVERLAP = 0.05;

    totalTrajectoryCount = 0;
    finalLoopQualities = [];
    loopCertainties = [];
    loopOverlap = [];
    for i = 1:max(loopCommunities)
        thisLoopIDs = find(loopCommunities == i);

        totalTrajectories = [];
        allPoints = [];
        allStates = [];

        notLoopIDs = 1:size(startConnectionsWeighted, 1);
        inLoopIDs = zeros(1, size(startConnectionsWeighted, 1));
        for j = 1:length(thisLoopIDs)
            loopSubtrajectoryID = thisLoopIDs(j);

            subLoops = allLoops{loopSubtrajectoryID};

            notLoopIDs(subLoops) = 0;
            inLoopIDs(subLoops) = 1;

            totalTrajectories(j) = sum(trajectoryCounts(subLoops));

            totalTrajectoryCount = totalTrajectoryCount + totalTrajectories(j);

            if ~USE_ALL_LOOPS
                allPoints = unique([allPoints trajectoryIndices{loopSubtrajectoryID}']);
                allStates = unique([allStates trajectoryStateMap{loopSubtrajectoryID}]);
            else
                for k = subLoops
                    allPoints = unique([allPoints trajectoryIndices{k}']);
                    allStates = unique([allStates trajectoryStateMap{k}]);
                end
            end
        end

        if USE_ALL_LOOPS
            inLoopIDs = find(inLoopIDs);% Use ALL loop trajectories in cluster
        else
            inLoopIDs = thisLoopIDs; % Use only trajectories in cluster
        end

        %%

        notLoopIDs(notLoopIDs == 0) = [];
        thisOverlaps = zeros(1, size(startConnectionsWeighted, 1));
        for j = 1:size(startConnectionsWeighted, 1)
            thisOverlaps(j) = sum(stateOverlap(j, notLoopIDs));
        end

        loopCertainties(i) = (sum(totalTrajectories));
        loopOverlap(i) = max(MIN_OVERLAP, sum(loopWeights(i,:) .* (thisOverlaps)));
        finalLoopQualities(i) = loopCertainties(i) / loopOverlap(i).^2;
    end
    
    thisLoopQuality = (sum(finalLoopQualities));
    
    loopTotalQualities(clusterCounts) = thisLoopQuality;
    
    if bestLoopQuality < thisLoopQuality
        bestLoopQuality = thisLoopQuality;
        
        bestClusterCount = clusterCounts;
    end

%     if length(loopQualities) > 1
%         [loopQualityClusters, centroids] = kmeans(log(loopQualities'), 2);
% 
%         bestThreshold = (exp(mean(centroids)));
%         [~, bestCentroid] = max(centroids);
%     else
%         loopQualityClusters = 1;
%         bestCentroid = 1;
%     end
    
end

[~, bestClusterCount] = max(diff(loopTotalQualities));
bestClusterCount = bestClusterCount + 1;

finalCommunities = cluster(loopLinkage,'maxclust',bestClusterCount);
goodLoops = 1:max(finalCommunities);

PLOT_LOOPS = goodLoops;

%%
% 
% % Calculate loop weights
% totalTrajectoryCount = 0;
% 
% loopWeights = [];
% for i = 1:max(finalCommunities)
%     thisLoopIDs = find(finalCommunities == i);
%     
%     totalTrajectories = [];
%     
%     notLoopIDs = 1:size(startConnectionsWeighted, 1);
%     inLoopIDs = zeros(1, size(startConnectionsWeighted, 1));
%     thisWeights = zeros(1, size(startConnectionsWeighted, 1));
%     for j = 1:length(thisLoopIDs)
%         loopSubtrajectoryID = thisLoopIDs(j);
%         
%         subLoops = allLoops{loopSubtrajectoryID};
%         
%         notLoopIDs(subLoops) = 0;
%         inLoopIDs(subLoops) = 1;
%         
%         thisWeights(subLoops) = thisWeights(subLoops) + 1;
%         
%         totalTrajectories(j) = sum(trajectoryCounts(subLoops));
%         
%         totalTrajectoryCount = totalTrajectoryCount + totalTrajectories(j);
%     end
%     thisWeights = thisWeights / sum(thisWeights);
%     loopWeights(i,:) = thisWeights;
%     
%     inLoopIDs = find(inLoopIDs);
% end
% 
% 
% % Get loop point weights
% loopPointWeights = zeros(length(PLOT_LOOPS), size(asymmetricProbabilities,1));
% for i = 1:max(loopCommunities)
%     subLoops = find(loopWeights(i,:) > 0);
%     
%     allStates = [];
%     for j = 1:length(subLoops)
%         thisTrajectory = subLoops(j);
%         thisStates = trajectoryStateMap{thisTrajectory};
%         
%         allStates = [allStates, thisStates];
%     end
%     
%     loopPointWeight = zeros(1, size(asymmetricProbabilities,1));
%     for j = 1:length(allStates)
%         trajectoryID = stateTrajectoryMap(allStates(j), 1);
%         
%         trajectoryWeight = loopWeights(i, trajectoryID);
%         
%         maxWeights = stateWeights(allStates(j),:);
%         maxWeights = maxWeights / max(maxWeights);
%         
%         loopPointWeight = loopPointWeight + maxWeights * trajectoryWeight;
%     end
%     
%     thisLoopPoints = find(loopPointWeight > 0.001);
%     loopPointWeights(i,thisLoopPoints) = loopPointWeight(thisLoopPoints);
%     
%     
%     figure(i);
%     clf;
%     hold on;
%     h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
%     h.Color(4) = 0.2;
%     %     scatter3(finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,1), finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,2), finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,3), 32, currentDistribution);
%         scatter3(finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,1), finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,2), finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,3), 32, loopPointWeight(thisLoopPoints));
%     colormap(parula(256));
%     subLoops = find(loopWeights(i,:) > 0);
%     
%     for j = subLoops
%         thisStates = trajectoryStateMap{j};
%     end
%     
% end
% 
% %%
% 
% figure(1);
% 
% for j = 57%1:length(loopOverlap)
% %     if loopCommunities(j) ~= 3
% %         continue;
% %     end
%     
%     clf;
%     hold on;
%     h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
%     h.Color(4) = 0.2;
%     subLoops = allLoops{j};
%     for i = 1:length(subLoops)    
%         plot3(trajectoryMeans{subLoops(i)}*pcaBasis(:,1), trajectoryMeans{subLoops(i)}*pcaBasis(:,2), trajectoryMeans{subLoops(i)}*pcaBasis(:,3), 'LineWidth', 3);
%     end
%     title(j);
%     colormap(parula(256));
% end
% 
% %% Assign poor loops to valid loops
% 
% stateOverlapNormalized = stateOverlap;
% loopOverlaps = [];
% for i = 1:length(loopQualityClusters)
%     thisLoopIDs = find(loopCommunities == i);
%     
%     for j = 1:length(loopQualityClusters)
%         maxOverlaps = [];
%         otherLoopIDs = find(loopCommunities == j);
%         
%         for k = 1:length(thisLoopIDs)
%             maxOverlaps(k) = max(stateOverlapNormalized(thisLoopIDs(k),otherLoopIDs));
%         end
%         
%         loopOverlaps(i,j) = mean(maxOverlaps);
%     end
% end
% loopOverlaps = loopOverlaps ./ max(loopOverlaps,[],2);
% 
% goodLoops = find(loopQualityClusters == bestCentroid);
% loopMerges = 1:length(loopQualityClusters);
% for i = 1:length(loopQualityClusters)
%     if loopQualityClusters(i) ~= bestCentroid
%         [~, bestMerge] = max(loopOverlaps(i,goodLoops));
%         
%         loopMerges(i) = goodLoops(bestMerge);
%     end
% end
% 
% newCommunities = loopCommunities;
% for i = 1:length(loopQualityClusters)
%     newCommunities(newCommunities == i) = loopMerges(i);
% end
% 
% loopRemap = unique(newCommunities);
% finalCommunities = newCommunities;
% for i = 1:length(loopRemap)
%     finalCommunities(newCommunities == loopRemap(i)) = i;
% end
% 
% PLOT_LOOPS = 1:max(finalCommunities);
% 
% goodLoops = 1:max(finalCommunities);
% 
% 


%% Recalculate loop weights
loopWeights = [];
for i = PLOT_LOOPS
    thisLoopIDs = find(finalCommunities == i);
    
    totalTrajectories = [];
    
    notLoopIDs = 1:size(startConnectionsWeighted, 1);
    inLoopIDs = zeros(1, size(startConnectionsWeighted, 1));
    thisWeights = zeros(1, size(startConnectionsWeighted, 1));
    for j = 1:length(thisLoopIDs)
        loopSubtrajectoryID = thisLoopIDs(j);
        
        subLoops = allLoops{loopSubtrajectoryID};
        
        notLoopIDs(subLoops) = 0;
        inLoopIDs(subLoops) = 1;
        
        thisWeights(subLoops) = thisWeights(subLoops) + 1;
        
        totalTrajectories(j) = sum(trajectoryCounts(subLoops));
        
        totalTrajectoryCount = totalTrajectoryCount + totalTrajectories(j);
    end
    thisWeights = thisWeights / sum(thisWeights);
    loopWeights(i,:) = thisWeights;
    
    inLoopIDs = find(inLoopIDs);
end



%% Get loop point weights
TEST_LOOP = 3;

loopPointWeights = zeros(length(PLOT_LOOPS), size(asymmetricProbabilities,1));
for i = PLOT_LOOPS
    subLoops = find(loopWeights(i,:) > 0);
    
    allStates = [];
    for j = 1:length(subLoops)
        thisTrajectory = subLoops(j);
        thisStates = trajectoryStateMap{thisTrajectory};
        
        allStates = [allStates, thisStates];
    end
    
    loopPointWeight = zeros(1, size(asymmetricProbabilities,1));
    for j = 1:length(allStates)
        trajectoryID = stateTrajectoryMap(allStates(j), 1);
        
        trajectoryWeight = loopWeights(i, trajectoryID);
        
        maxWeights = stateWeights(allStates(j),:);
        maxWeights = maxWeights / max(maxWeights);
        
        loopPointWeight = loopPointWeight + maxWeights * trajectoryWeight;
    end
    
    thisLoopPoints = find(loopPointWeight > 0.001);
    loopPointWeights(i,thisLoopPoints) = loopPointWeight(thisLoopPoints);
%     [~, startPoint] = max(thisLoopWeights{i});
%     %     startPoint = 1;
%     bestPoint = startPoint;
    
    
    
    %     lastValue = 0;
    %     rotateCount = 0;
    %     usedPoints = [];
    %     meanPoints = [];
    %     for j = 1:400
    %         currentDistribution = asymmetricProbabilities(thisLoopPoints(bestPoint), thisLoopPoints);
    % %         currentDistribution = sum(currentDistribution * asymmetricProbabilities(thisLoopPoints, thisLoopPoints)) .* thisLoopWeights{i};
    %         currentDistribution = currentDistribution / sum(currentDistribution);
    %
    %         if ~ismember(bestPoint, usedPoints)
    %             meanPoints(1+size(meanPoints,1),:) = currentDistribution * finalDynamicsStream(thisLoopPoints,:);
    %             lastValue = bestPoint;
    %             rotateCount = 0;
    %             usedPoints = [bestPoint usePoints];
    %         else
    %             rotateCount = rotateCount + 1;
    %         end
    % %         currentPoints = thisLoopPoints(find(currentDistribution > 0.001));
    % %         currentPoints = currentPoints + 1;
    % %         currentPoints = intersect(currentPoints, thisLoopPoints);
    %
    %         newDistribution = [zeros(1,rotateCount) currentDistribution];
    %         newDistribution(end-rotateCount+1:end) = [];
    %         currentDistribution = newDistribution;
    %         currentDistribution = currentDistribution / sum(currentDistribution);
    %
    %         currentDistribution = currentDistribution * asymmetricProbabilities(thisLoopPoints, thisLoopPoints);
    %         currentDistribution = currentDistribution / sum(currentDistribution);
    %
    %         bestOverLap = currentDistribution * asymmetricProbabilities(thisLoopPoints, thisLoopPoints);
    %         [~, bestPoint] = max(bestOverLap);
    %     end
    
    %     bestPoint = startPoint;
    %
    %     lastValue = 0;
    %     rotateCount = 0;
    %     backwardMeanPoints = [];
    %     for j = 1:200
    %         currentDistribution = asymmetricProbabilities(thisLoopPoints, thisLoopPoints(bestPoint))';
    % %         currentDistribution = sum(currentDistribution * asymmetricProbabilities(thisLoopPoints, thisLoopPoints)) .* thisLoopWeights{i};
    %         currentDistribution = currentDistribution / sum(currentDistribution);
    %
    %         if lastValue ~= bestPoint
    %             backwardMeanPoints(1+size(backwardMeanPoints,1),:) = currentDistribution * finalDynamicsStream(thisLoopPoints,:);
    %             lastValue = bestPoint;
    %             rotateCount = 0;
    %         else
    %             rotateCount = rotateCount + 1;
    %         end
    % %         currentPoints = thisLoopPoints(find(currentDistribution > 0.001));
    % %         currentPoints = currentPoints + 1;
    % %         currentPoints = intersect(currentPoints, thisLoopPoints);
    %
    %         newDistribution = [currentDistribution zeros(1,rotateCount) ];
    %         newDistribution(1:rotateCount) = [];
    %         currentDistribution = newDistribution;
    %         currentDistribution = currentDistribution / sum(currentDistribution);
    %
    %         currentDistribution = currentDistribution * asymmetricProbabilities(thisLoopPoints, thisLoopPoints)';
    %         currentDistribution = currentDistribution / sum(currentDistribution);
    %
    %         bestOverLap = currentDistribution * asymmetricProbabilities(thisLoopPoints, thisLoopPoints)';
    %         [~, bestPoint] = max(bestOverLap);
    %     end
    
    %     bestPoint = startPoint;
    %
    %     backwardMeanPoints = [];
    %     for j = 1:20
    %         currentDistribution = asymmetricProbabilities(thisLoopPoints, thisLoopPoints(bestPoint))';
    % %         currentDistribution = sum(currentDistribution * asymmetricProbabilities(thisLoopPoints, thisLoopPoints)) .* thisLoopWeights{i};
    %         currentDistribution = currentDistribution / sum(currentDistribution);
    %
    %         backwardMeanPoints(j,:) = currentDistribution * finalDynamicsStream(thisLoopPoints,:);
    %
    % %         currentPoints = thisLoopPoints(find(currentDistribution > 0.001));
    % %         currentPoints = currentPoints + 1;
    % %         currentPoints = intersect(currentPoints, thisLoopPoints);
    %
    %         lerpValue = 0;
    %
    %         newDistribution = [currentDistribution 0 ];
    %         newDistribution(1) = [];
    %         currentDistribution = currentDistribution * (1-lerpValue) + newDistribution * lerpValue;
    %         currentDistribution = currentDistribution / sum(currentDistribution);
    %
    %         currentDistribution = currentDistribution * asymmetricProbabilities(thisLoopPoints, thisLoopPoints)';
    %         currentDistribution = currentDistribution / sum(currentDistribution);
    %
    %         bestOverLap = currentDistribution * asymmetricProbabilities(thisLoopPoints, thisLoopPoints)';
    %         [~, bestPoint] = max(bestOverLap);
    %     end
    
    figure(i);
    clf;
    hold on;
    h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
    h.Color(4) = 0.2;
    %     scatter3(finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,1), finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,2), finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,3), 32, currentDistribution);
        scatter3(finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,1), finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,2), finalDynamicsStream(thisLoopPoints,:)*pcaBasis(:,3), 32, loopPointWeight(thisLoopPoints));
    colormap(parula(256));
    %     plot3(meanPoints*pcaBasis(:,1),meanPoints*pcaBasis(:,2),meanPoints*pcaBasis(:,3), 'LineWidth', 2);
    %     plot3(backwardMeanPoints*pcaBasis(:,1),backwardMeanPoints*pcaBasis(:,2),backwardMeanPoints*pcaBasis(:,3), 'LineWidth', 2);
    subLoops = find(loopWeights(i,:) > 0);
    
    for j = subLoops
        thisStates = trajectoryStateMap{j};
        
        %                 scatter3(trajectoryMeans(:,:,j)*pcaBasis(:,1), trajectoryMeans(:,:,j)*pcaBasis(:,2), trajectoryMeans(:,:,j)*pcaBasis(:,3), 32, ones(size(thisStates)) * loopOrdering(j));
%         scatter3(trajectoryMeans{j}*pcaBasis(:,1), trajectoryMeans{j}*pcaBasis(:,2), trajectoryMeans{j}*pcaBasis(:,3), 32, loopStateOrders(i,thisStates));
        %         scatter3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 32, loopPointOrders(i,:));
    end
    
end



%%

NEAREST_ANGLES = 5;
USE_ALL_LOOPS = 1;
MIN_OVERLAP = 0.05;
USE_FIXED_POINTS = 0;
if USE_FIXED_POINTS
    USE_ALL_LOOPS = 0;
end
TEMPORAL_NEIGHBORHOODS = [5, 5, 10, 5, 5, 5, 5];
PHASE_BINS = 50;
PHASE_THRESHOLD = 1/20;
FORCE_LOOP = 1;
MIN_WEIGHT = exp(-2);
PHASE_MIN_WEIGHT = exp(-1);
    
useLoopPointWeights = loopPointWeights;

rawDistances = asymmetricProbabilities ./ max(asymmetricProbabilities);

totalTrajectoryCount = 0;
loopMeans = {};
loopSTDs = {};
loopVelocities = {};
loopSpeeds = {};
loopPointWeightings = [];
loopIndicies = {};
loopPhases = [];
for i = PLOT_LOOPS
    [~, bestTrajectory] = max(loopWeights(i,:));
    
    thisLoopIDs = find(finalCommunities == i);
    %             thisLoopIDs(1) = [];
    %             thisLoopIDs(end) = [];
    
    if isempty(thisLoopIDs)
        continue;
    end
    
    totalTrajectories = [];
    allPoints = [];
    allStates = [];
    
    notLoopIDs = 1:size(startConnectionsWeighted, 1);
    inLoopIDs = zeros(1, size(startConnectionsWeighted, 1));
    for j = 1:length(thisLoopIDs)
        loopSubtrajectoryID = thisLoopIDs(j);
        
        subLoops = allLoops{loopSubtrajectoryID};
        
        notLoopIDs(subLoops) = 0;
        inLoopIDs(subLoops) = 1;
        
        totalTrajectories(j) = sum(trajectoryCounts(subLoops));
        
        totalTrajectoryCount = totalTrajectoryCount + totalTrajectories(j);
        
        if ~USE_ALL_LOOPS
            allPoints = unique([allPoints trajectoryIndices{loopSubtrajectoryID}']);
%             allStates = unique([allStates trajectoryStateMap{loopSubtrajectoryID}]);
        else
            for k = subLoops
                allPoints = unique([allPoints trajectoryIndices{k}']);
%                 allStates = unique([allStates trajectoryStateMap{k}]);
            end
        end
    end
    
    if USE_ALL_LOOPS
        inLoopIDs = find(inLoopIDs);% Use ALL loop trajectories in cluster
    else
        inLoopIDs = thisLoopIDs; % Use only trajectories in cluster
    end
    
    useLoopPointWeights(i, setdiff(1:size(useLoopPointWeights,2), allPoints)) = 0;
    allPoints = intersect(allPoints, find(useLoopPointWeights(i,:) > MIN_WEIGHT));
    
    %% Get local phase
    
    averagePoints = allPoints;
    averagePoints = averagePoints' + (-TAIL_LENGTH:TAIL_LENGTH);
    averagePoints(averagePoints < 1) = [];
    averagePoints(averagePoints > size(asymmetricProbabilities,1)) = [];
    averagePoints = unique(averagePoints)';
    
    loopProbabilities = asymmetricProbabilities(averagePoints,averagePoints);
    loopProbabilities = loopProbabilities ./ sum(loopProbabilities,2);
    
    % Get detailed balance decomp
    [steadyState, ~] = eigs(loopProbabilities', 1);
    fluxMatrix = diag(steadyState) * loopProbabilities;
    symmetricMatrix = (fluxMatrix+fluxMatrix.')/2;
    antisymmetricMatrix = (fluxMatrix-fluxMatrix.')/2;
    symmetricMatrix = symmetricMatrix ./ steadyState;
    antisymmetricMatrix = antisymmetricMatrix ./ steadyState;
    
    % Rescale diffusion map
    %         diffusedProbabilities = max(0, expm(symmetricMatrix) - eye(size(symmetricMatrix)) + antisymmetricMatrix);
    diffusedProbabilities = max(0, symmetricMatrix + sign(antisymmetricMatrix).*abs(antisymmetricMatrix).^(1/8));
    %     diffusedProbabilities = antisymmetricMatrix;
    %     diffusedProbabilities = max(0, symmetricMatrix + antisymmetricMatrix);
    % diffusedProbabilities = symmetricMatrix + antisymmetricMatrix;
    diffusedProbabilities = diffusedProbabilities ./ sum(diffusedProbabilities,2);
    
    [eigenVectors, eigenValues] = eigs(diffusedProbabilities, min(10, length(averagePoints)));
    eigenValues = diag(eigenValues);
    
    [~, bestComplexEigenVector] = max((abs(imag(eigenValues)) > 0) .* abs(eigenValues));
    
    loopPointPhase = angle(eigenVectors(:,bestComplexEigenVector));
    
    inLoopIndices = find(ismember(averagePoints, allPoints));
    
%     loopPhases{i} = loopPointPhase(inLoopIndices);
    loopIndicies{i} = allPoints;
    
    figure(5);
    clf;
    hold on;
    h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
    h.Color(4) = 0.2;
    colormap(hsv(256));
    for j = PLOT_LOOPS
        plotIndices = averagePoints(inLoopIndices);
        loopSpaceIndicies = find(ismember(averagePoints, plotIndices));

        scatter3(finalDynamicsStream(plotIndices,:)*pcaBasis(:,1), finalDynamicsStream(plotIndices,:)*pcaBasis(:,2), finalDynamicsStream(plotIndices,:)*pcaBasis(:,3), 32, loopPointPhase(loopSpaceIndicies));
    end
    %             caxis([-pi pi]);
    colorbar
    
    %% Unwrap temporally
    
    startTrajectory = thisLoopIDs(1);
    %             if numFixedPoints == 1
    %                 startTrajectory = bestTrajectory;
    %             end
    
    thisIndicies = trajectoryIndices{startTrajectory} + 1;
    thisIndicies = setdiff(trajectoryIndices{startTrajectory}, thisIndicies);
    thisIndicies(thisIndicies > size(dataTrajectoryMap,2)) = [];
    phaseStarts = thisIndicies;
    
    thisIndicies = inLoopIndices + 1;
    thisIndicies = setdiff(inLoopIndices, thisIndicies);
    thisIndicies(thisIndicies > size(dataTrajectoryMap,2)) = [];
    loopStartPoints = averagePoints(thisIndicies);
    if size(loopStartPoints,2) > size(loopStartPoints,1)
        loopStartPoints = loopStartPoints';
    end
    
    thisIndicies = inLoopIndices - 1;
    thisIndicies = setdiff(inLoopIndices, thisIndicies);
    thisIndicies(thisIndicies > size(dataTrajectoryMap,2)) = [];
    loopEndPoints = averagePoints(thisIndicies);
    if size(loopEndPoints,2) > size(loopEndPoints,1)
        loopEndPoints = loopEndPoints';
    end
    for j = 1:length(phaseStarts)
        if ~ismember(phaseStarts(j), loopStartPoints)
            loopStartPoints = [loopStartPoints; phaseStarts(j)];
            loopEndPoints = [loopEndPoints; phaseStarts(j)-1];
        end
    end
    
    discontinuities = intersect(find(~stateHasNext), averagePoints(inLoopIndices));
    for j = 1:length(discontinuities)
        if ~ismember(discontinuities(j), loopEndPoints)
            loopStartPoints = [loopStartPoints; discontinuities(j)+1];
            loopEndPoints = [loopEndPoints; discontinuities(j)];
        end
    end
    
    loopStartPoints = sort(loopStartPoints);
    loopEndPoints = sort(loopEndPoints);
    
    allLoopStarts = [];
    for j = 1:length(phaseStarts)
        thisTraceIndex = find(phaseStarts(j) == loopStartPoints);
        thisTracesIDs = loopStartPoints(thisTraceIndex);
        thisTracesIDsLoopSpace = find(ismember(averagePoints, thisTracesIDs));
        
        allLoopStarts = [allLoopStarts, thisTracesIDsLoopSpace];
    end
    
    averageStartPhase = angle(mean(exp(1i*loopPointPhase(allLoopStarts))));
    loopPointPhase = wrapToPi(loopPointPhase - averageStartPhase);
    
    unwrappedPhases = nan(size(asymmetricProbabilities,1), 1);
    exploredPoints = [];
    for j = 1:length(phaseStarts)
        thisTraceIndex = find(phaseStarts(j) == loopStartPoints);
        thisTracesIDs = loopStartPoints(thisTraceIndex)-(0:TEMPORAL_NEIGHBORHOODS(i));%:loopEndPoints(thisTraceIndex);
        thisTracesIDs(thisTracesIDs > size(asymmetricProbabilities,1)) = [];
        thisTracesIDs(thisTracesIDs < 1) = [];
        thisTracesIDsLoopSpace = find(ismember(averagePoints, thisTracesIDs));
        
        unwrappedPhases(thisTracesIDs) = 0:0+length(thisTracesIDs)-1;
        exploredPoints = [exploredPoints, thisTracesIDs];
        
%         figure(5);
%         clf;
%         hold on;
%         h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
%         h.Color(4) = 0.2;
%         colormap(parula(256));
%         for j = PLOT_LOOPS
%             plotIndices = averagePoints(inLoopIndices);
%             loopSpaceIndicies = find(ismember(averagePoints, plotIndices));
%             
%             scatter3(finalDynamicsStream(plotIndices,:)*pcaBasis(:,1), finalDynamicsStream(plotIndices,:)*pcaBasis(:,2), finalDynamicsStream(plotIndices,:)*pcaBasis(:,3), 32, unwrappedPhases(plotIndices));
%         end
%         %             caxis([-pi pi]);
%         colorbar
    end
    
    pointsToMerge = setdiff(averagePoints(inLoopIndices), exploredPoints);
    waitHandle = parfor_progressbar(length(pointsToMerge) + length(averagePoints), ['Unwrapping phase (' num2str(i) ' of ' num2str(length(PLOT_LOOPS)) ')']);
    while length(pointsToMerge) > 0
        hasPhase = zeros(size(unwrappedPhases));
        hasPhase(~isnan(unwrappedPhases)) = 1;
        bestMerge = asymmetricProbabilities(pointsToMerge,:) * hasPhase .* useLoopPointWeights(i,pointsToMerge)';
        %                 bestMerge = rawDistances(pointsToMerge,:) * hasPhase;
        [~, bestMerge] = max(bestMerge);
        
        thisTracesIDs = pointsToMerge(bestMerge)-(0:TEMPORAL_NEIGHBORHOODS(i));%:loopEndPoints(thisTraceIndex);
        thisTracesIDs(thisTracesIDs > size(asymmetricProbabilities,1)) = [];
        thisTracesIDs(thisTracesIDs < 1) = [];
        thisTracesIDsLoopSpace = find(ismember(averagePoints, thisTracesIDs));
        
        originalPhase = loopPointPhase(thisTracesIDsLoopSpace(1));
        
        thisWeights = asymmetricProbabilities(pointsToMerge(bestMerge),:);% .* useLoopPointWeights(i,:);
        validWeights = find(~isnan(thisWeights .* unwrappedPhases') & thisWeights > 0);
        [sortedWeights, sortOrder] = sort(thisWeights(validWeights), 'descend');
        topWeights = sortOrder;%(1:NEAREST_ANGLES);
        useIDs = validWeights(topWeights);
        useWeights = thisWeights(useIDs);
        useWeights = useWeights / sum(useWeights);
        bestStartPhase = nansum(useWeights .* unwrappedPhases(useIDs)');
        
        if sum(isnan(unwrappedPhases(thisTracesIDs))) < length(thisTracesIDs)
            realValues = find(~isnan(unwrappedPhases(thisTracesIDs)));
            [addedValue, addedIndex] = max(unwrappedPhases(thisTracesIDs(realValues)));
            fitPhase = fit([thisTracesIDs(realValues), thisTracesIDs(realValues(addedIndex))-1]', [unwrappedPhases(thisTracesIDs(realValues))',addedValue+1]', 'poly1');
            
            nanValues = find(isnan(unwrappedPhases(thisTracesIDs)));
            unwrappedPhases(thisTracesIDs(nanValues)) = feval(fitPhase, thisTracesIDs(nanValues));
        else
            if length(thisTracesIDs) == 1
                unwrappedPhases(thisTracesIDs) = bestStartPhase;
            else
                unwrappedPhases(thisTracesIDs) = bestStartPhase:bestStartPhase+length(thisTracesIDs)-1;
            end
        end
        
        
        removePoints = [thisTracesIDs];
        
        pointsToMerge(ismember(pointsToMerge, removePoints)) = [];
        
        
%         if abs(phaseOffset) > 0
%             test=1;
%         end
%         
%         if max(unwrap(loopPointPhase(thisTracesIDsLoopSpace)) + phaseOffset) > pi
%             test=1;
%         end

        waitHandle.iterate(1);
    end
    
    % Fix original points
    [angleCounts, edges] = histcounts(unwrappedPhases(averagePoints(inLoopIndices)), PHASE_BINS);
    maxTrajectorires = max(trajectoryCounts(thisLoopIDs));
    
    angleWeights = [];
    for j = 1:length(angleCounts)
        thisAngles = find(unwrappedPhases(averagePoints(inLoopIndices)) >= edges(j) & unwrappedPhases(averagePoints(inLoopIndices)) < edges(j+1));
        
        angleWeights(j) = sum(useLoopPointWeights(i,thisAngles));
    end
    angleWeights = filterData(angleWeights, 2);
    
    badPhases = unique([1, find(angleWeights < maxTrajectorires * PHASE_THRESHOLD), length(angleWeights)]);
    splitPoint = find(diff(badPhases) > 2);
    splitPoint = splitPoint(1);
    
%     minPhase = edges(badPhases(splitPoint) + 1);
%     maxPhase = edges(badPhases(splitPoint + 1));
%     minPhase = min(unwrappedPhases(averagePoints(inLoopIndices)));%quantile(unwrappedPhases, 0.01);
%     maxPhase = max(unwrappedPhases(averagePoints(inLoopIndices)));%quantile(unwrappedPhases, 0.99);

    goodPoints = find(useLoopPointWeights(i,averagePoints(inLoopIndices)) > PHASE_MIN_WEIGHT);

    minPhase = quantile(unwrappedPhases(averagePoints(inLoopIndices(goodPoints))), 0.01);
    maxPhase = quantile(unwrappedPhases(averagePoints(inLoopIndices(goodPoints))), 0.99);
    
    if ~FORCE_LOOP
        unwrappedPhases = (unwrappedPhases - minPhase) / (maxPhase - minPhase) * 2*pi - pi;
    else
        unwrappedPhases = wrapToPi((unwrappedPhases - minPhase) / (maxPhase - minPhase) * 2*pi - pi);
    end
    
    pointsToMerge = averagePoints;%(inLoopIndices);
    while length(pointsToMerge) > 0
        hasPhase = ones(size(unwrappedPhases));
%         hasPhase(~isnan(unwrappedPhases)) = 1;
%         hasPhase(exploredPoints) = 0;
        bestMerge = asymmetricProbabilities(pointsToMerge,:) * hasPhase .* useLoopPointWeights(i,pointsToMerge)';
        %                 bestMerge = rawDistances(pointsToMerge,:) * hasPhase;
        [~, bestMerge] = max(bestMerge);
        
        thisTracesIDs = pointsToMerge(bestMerge);
        thisTracesIDsLoopSpace = find(ismember(averagePoints, thisTracesIDs));
        
        originalPhase = loopPointPhase(thisTracesIDsLoopSpace(1));
        
        thisWeights = asymmetricProbabilities(pointsToMerge(bestMerge),:);% .* useLoopPointWeights(i,:);
        validWeights = find(~isnan(thisWeights .* unwrappedPhases') & thisWeights > 0);
        [sortedWeights, sortOrder] = sort(thisWeights(validWeights), 'descend');
        topWeights = sortOrder;%(1:NEAREST_ANGLES);
        useIDs = validWeights(topWeights);
        useWeights = thisWeights(useIDs);
        useWeights = useWeights / sum(useWeights);
        if FORCE_LOOP
            bestStartPhase = angle(nansum(useWeights .* exp(1i*unwrappedPhases(useIDs)')));
        else
            bestStartPhase = nansum(useWeights .* unwrappedPhases(useIDs)');
        end
        
        
        unwrappedPhases(thisTracesIDs) = bestStartPhase;
        
        removePoints = [thisTracesIDs];
        
        pointsToMerge(ismember(pointsToMerge, removePoints)) = [];
        %
        %                 figure(5+currentFixedPoint);
        %                 clf;
        %                 hold on;
        %                 h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
        %                 h.Color(4) = 0.2;
        %                 colormap(parula(256));
        %                 for j = PLOT_LOOPS
        %                     plotIndices = averagePoints(inLoopIndices);
        %                     loopSpaceIndicies = find(ismember(averagePoints, plotIndices));
        %
        %                     scatter3(finalDynamicsStream(plotIndices,:)*pcaBasis(:,1), finalDynamicsStream(plotIndices,:)*pcaBasis(:,2), finalDynamicsStream(plotIndices,:)*pcaBasis(:,3), 32, unwrappedPhases(plotIndices));
        %                 end
        
%         if abs(phaseOffset) > 0
%             test=1;
%         end
%         
%         if max(unwrap(loopPointPhase(thisTracesIDsLoopSpace)) + phaseOffset) > pi
%             test=1;
%         end
        waitHandle.iterate(1);
    end
    close(waitHandle);
    % Finalize phase
%     [angleCounts, edges] = histcounts(unwrappedPhases(averagePoints(inLoopIndices)), PHASE_BINS);
%     maxTrajectorires = max(trajectoryCounts(thisLoopIDs));
%     
%     angleWeights = [];
%     for j = 1:length(angleCounts)
%         thisAngles = find(unwrappedPhases(averagePoints(inLoopIndices)) >= edges(j) & unwrappedPhases(averagePoints(inLoopIndices)) < edges(j+1));
%         
%         angleWeights(j) = sum(useLoopPointWeights(i,thisAngles));
%     end
%     angleWeights = filterData(angleWeights, 2);
%     
%     badPhases = unique([1, find(angleWeights < maxTrajectorires * PHASE_THRESHOLD), length(angleWeights)]);
%     splitPoint = find(diff(badPhases) > 2);
%     splitPoint = splitPoint(1);
%     
%     minPhase = edges(badPhases(splitPoint) + 1);
%     maxPhase = edges(badPhases(splitPoint + 1));
%     minPhase = min(unwrappedPhases(averagePoints(inLoopIndices)));%quantile(unwrappedPhases, 0.01);
%     maxPhase = max(unwrappedPhases(averagePoints(inLoopIndices)));%quantile(unwrappedPhases, 0.99);
    minPhase = quantile(unwrappedPhases(averagePoints(inLoopIndices)), 0.01);
    maxPhase = quantile(unwrappedPhases(averagePoints(inLoopIndices)), 0.99);
    
    if ~FORCE_LOOP
        unwrappedPhases = (unwrappedPhases - minPhase) / (maxPhase - minPhase) * 2*pi - pi;
    else
        unwrappedPhases = wrapToPi((unwrappedPhases - minPhase) / (maxPhase - minPhase) * 2*pi - pi);
    end
    
    loopPhases(i,:) = unwrappedPhases;
    
    %             figure(6);
    %             clf;
    %             hold on;
    %             h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
    %             h.Color(4) = 0.2;
    %             colormap(hsv(256));
    %             for j = PLOT_LOOPS
    %                 plotIndices = averagePoints(inLoopIndices);
    %                 loopSpaceIndicies = find(ismember(averagePoints, plotIndices));
    %
    %                 scatter3(finalDynamicsStream(plotIndices,:)*pcaBasis(:,1), finalDynamicsStream(plotIndices,:)*pcaBasis(:,2), finalDynamicsStream(plotIndices,:)*pcaBasis(:,3), 32, loopPointPhase(loopSpaceIndicies));
    %             end
    % %             caxis([-pi pi]);
    %             colorbar
    
    figure(5+i);
    clf;
    hold on;
    h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
    h.Color(4) = 0.2;
    colormap(parula(256));
    for j = PLOT_LOOPS
        plotIndices = averagePoints(inLoopIndices);
        loopSpaceIndicies = find(ismember(averagePoints, plotIndices));
        
        scatter3(finalDynamicsStream(plotIndices,:)*pcaBasis(:,1), finalDynamicsStream(plotIndices,:)*pcaBasis(:,2), finalDynamicsStream(plotIndices,:)*pcaBasis(:,3), 32, unwrappedPhases(plotIndices));
    end
    %             caxis([-pi pi]);
    colormap(hsv(256));
    colorbar
    drawnow();
end

%% Approximate loops


NUM_SEGMENTS = 100;
loopPhase = -pi:2*pi/NUM_SEGMENTS:pi;

loopMeans = {};
loopSTDs = {};
loopPointWeightings = [];
loopVelocities = {};
loopSpeeds = {};
for i = PLOT_LOOPS    
    
%     [phaseValues, sortedPhases] = sort(loopPhases);
%     newStart = find(phaseValues == startPhase);
%     sortedPhases = [sortedPhases(newStart:end) sortedPhases(1:newStart-1)];
%     
    loopPhase = loopPhase;
    splinePhases = loopPhases(i,:);
%     splineWeights = loopPointWeight';
    splineWeights = useLoopPointWeights(i,:)'.^2;
    splineWeights = splineWeights / sum(splineWeights);
%     splineWeights = sqrt(splineWeights);
%     splinePositions = finalDynamicsStream;
    
    LOOP_SIGMA = pi/20;
    
    loopMean = [];
    loopSTD = [];
    loopWeight = [];
    loopVelocity = [];
    loopSpeed = [];
    for j = 1:length(loopPhase)
        weights = exp(-(angleDiff(splinePhases, loopPhase(j))).^2/(2*LOOP_SIGMA^2))';
%         weights = exp(-(splinePhases - loopPhase(j)).^2/(2*LOOP_SIGMA^2));
        weights = weights / nansum(weights);
        
        loopWeight(j,:) = weights;
        
        weights = weights .* splineWeights;
        weights = weights / nansum(weights);
        
        nanValues = find(~isnan(weights));
        
        loopMean(j,:) = nansum(weights(nanValues) .* finalDynamicsStream(nanValues,:));
        loopSTD(j,:) = sqrt(var(finalDynamicsStream(nanValues,:), weights(nanValues)));
    end
    loopVelocity = [diff(loopMean, 1); zeros(1, size(loopMean,2))];
    loopSpeed = sqrt(sum(loopVelocity.^2,2));
    
    loopMeans{i} = loopMean;
    loopSTDs{i} = loopMean;
    loopPointWeightings(i,:,:) = loopWeight;
    loopVelocities{i} = loopVelocity;
    loopSpeeds{i} = loopSpeed;
    
    %%
    
    
    
    
    %%
    
    figure(10+i);
    clf;
    hold on;
    h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
    h.Color(4) = 0.2;
    colormap(parula(256));
    %             for j = PLOT_LOOPS
    %                 plotPath = [loopSplines{i}];
    subLoops = find(loopWeights(i,:) > 0);
    
    %                     for j = subLoops
    %                 %         h = plot3(trajectoryMeans(:,:,j)*pcaBasis(:,1), trajectoryMeans(:,:,j)*pcaBasis(:,2), trajectoryMeans(:,:,j)*pcaBasis(:,3), 'k', 'LineWidth', 3);
    %                 %         h.Color(4) = loopWeights(i,j) / max(loopWeights(i,:));
    %                         thisStates = trajectoryStateMap(j,:);
    %
    %                         scatter3(trajectoryMeans(:,:,j)*pcaBasis(:,1), trajectoryMeans(:,:,j)*pcaBasis(:,2), trajectoryMeans(:,:,j)*pcaBasis(:,3), 32, loopStateOrders(i,thisStates));
    %
    %                     end
    
    %                 scatter3(splinePositions*pcaBasis(:,1), splinePositions*pcaBasis(:,2), splinePositions*pcaBasis(:,3), 32, log(splineWeights));
    %                 scatter3(splinePositions*pcaBasis(:,1), splinePositions*pcaBasis(:,2), splinePositions*pcaBasis(:,3), 32, splinePhases);
    %                 scatter3(loopSplines{i}'*pcaBasis(:,1), loopSplines{i}'*pcaBasis(:,2), loopSplines{i}'*pcaBasis(:,3), 32, loopSplinePhases{i});
    
%     scatter3(finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,1), finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,2), finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,3), 32, splineWeights(loopIndicies{i}));
%     scatter3(finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,1), finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,2), finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,3), 32, unwrappedPhases(loopIndicies{i}));
    %                     scatter3(finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,1), finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,2), finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,3), 32, loopSplineWeights{i});
%     plot3(loopTrajectoryMeans*pcaBasis(:,1), loopTrajectoryMeans*pcaBasis(:,2), loopTrajectoryMeans*pcaBasis(:,3), 'Color', colors(colorValues(i),:), 'LineWidth', 3);
%         plot3(loopSplines{i}'*pcaBasis(:,1), loopSplines{i}'*pcaBasis(:,2), loopSplines{i}'*pcaBasis(:,3), 'LineWidth', 3);
        scatter3(loopMeans{i}*pcaBasis(:,1), loopMeans{i}*pcaBasis(:,2), loopMeans{i}*pcaBasis(:,3), 32, loopSpeeds{i});
    %             end
    
    drawnow();
end

%         figure(5);
%         clf;
%         hold on;
%         h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
%         h.Color(4) = 0.2;
%         for j = PLOT_LOOPS
%             plotPath = [loopSplines{i} loopSplines{i}(:,1)];
%             subLoops = find(loopWeights(i,:) > 0);
%
%             %     for j = subLoops
%             % %         h = plot3(trajectoryMeans(:,:,j)*pcaBasis(:,1), trajectoryMeans(:,:,j)*pcaBasis(:,2), trajectoryMeans(:,:,j)*pcaBasis(:,3), 'k', 'LineWidth', 3);
%             % %         h.Color(4) = loopWeights(i,j) / max(loopWeights(i,:));
%             %         thisStates = trajectoryStateMap(j,:);
%             %
%             %         scatter3(trajectoryMeans(:,:,j)*pcaBasis(:,1), trajectoryMeans(:,:,j)*pcaBasis(:,2), trajectoryMeans(:,:,j)*pcaBasis(:,3), 32, loopStateOrders(i,thisStates));
%             %
%             %     end
%
%             %     scatter3(splinePositions*pcaBasis(:,1), splinePositions*pcaBasis(:,2), splinePositions*pcaBasis(:,3), 32, loopStatePhases{i});
%             scatter3(loopSplines{i}'*pcaBasis(:,1), loopSplines{i}'*pcaBasis(:,2), loopSplines{i}'*pcaBasis(:,3), 32, loopSplinePhases{i});
%
%             scatter3(finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,1), finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,2), finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,3), 32, loopPhases{i});
%             %     scatter3(finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,1), finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,2), finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,3), 32, loopSplineWeights{i});
%             %     plot3(plotPath'*pcaBasis(:,1), plotPath'*pcaBasis(:,2), plotPath'*pcaBasis(:,3), 'Color', colors(colorValues(i),:), 'LineWidth', 3);
%             %     scatter3(loopSplines{i}*pcaBasis(:,1), loopSplines{i}*pcaBasis(:,2), loopSplines{i}*pcaBasis(:,3), 32, 'MarkerEdgeColor', colors(colorValues(i),:));
%         end

% %%

% loopQualities = loopQualities / sum(loopCertainties) * length(loopQualities);
% loopQualities(loopQualities > 20) = 20;

% if length(loopQualities) > 1
%     [loopQualityClusters, centroids] = kmeans(log(loopQualities'), 2);
%
%     bestThreshold = (exp(mean(centroids)));
%     [~, bestCentroid] = max(centroids);
% else
%     loopQualityClusters = 1;
%     bestCentroid = 1;
% end


% figure(5);
% clf;
% hold on;
% h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
% h.Color(4) = 0.2;
% for i = PLOT_LOOPS
%     plotPath = [loopSplines{i} loopSplines{i}(:,1)];
%     subLoops = find(loopWeights(i,:) > 0);
%
%     %     for j = subLoops
%     % %         h = plot3(trajectoryMeans(:,:,j)*pcaBasis(:,1), trajectoryMeans(:,:,j)*pcaBasis(:,2), trajectoryMeans(:,:,j)*pcaBasis(:,3), 'k', 'LineWidth', 3);
%     % %         h.Color(4) = loopWeights(i,j) / max(loopWeights(i,:));
%     %         thisStates = trajectoryStateMap(j,:);
%     %
%     %         scatter3(trajectoryMeans(:,:,j)*pcaBasis(:,1), trajectoryMeans(:,:,j)*pcaBasis(:,2), trajectoryMeans(:,:,j)*pcaBasis(:,3), 32, loopStateOrders(i,thisStates));
%     %
%     %     end
%
%     %     scatter3(splinePositions*pcaBasis(:,1), splinePositions*pcaBasis(:,2), splinePositions*pcaBasis(:,3), 32, loopStatePhases{i});
%     scatter3(loopSplines{i}'*pcaBasis(:,1), loopSplines{i}'*pcaBasis(:,2), loopSplines{i}'*pcaBasis(:,3), 32, loopSplinePhases{i});
%
%     scatter3(finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,1), finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,2), finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,3), 32, loopPhases{i});
%     %     scatter3(finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,1), finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,2), finalDynamicsStream(loopIndicies{i},:)*pcaBasis(:,3), 32, loopSplineWeights{i});
%     %     plot3(plotPath'*pcaBasis(:,1), plotPath'*pcaBasis(:,2), plotPath'*pcaBasis(:,3), 'Color', colors(colorValues(i),:), 'LineWidth', 3);
%     %     scatter3(loopSplines{i}*pcaBasis(:,1), loopSplines{i}*pcaBasis(:,2), loopSplines{i}*pcaBasis(:,3), 32, 'MarkerEdgeColor', colors(colorValues(i),:));
% end

%% Find emprical transitions
SIMILARITY_CUTOFF = 0.7;

neighborhoodQuality = inCounts ./ outCounts;

validWeights = loopPointWeightings(goodLoops,:,:);
validWeights = validWeights(:,1:end-1,:);
% linearWeights = reshape(validWeights, [size(validWeights,1)*size(validWeights,2), size(validWeights,3)]);
% linearWeights = linearWeights ./ sum(linearWeights,2);

% loopPointWeightsNormalized = loopPointWeights ./ sqrt(sum(loopPointWeights.^2,2));
% loopPointWeightsNormalized = loopPointWeights ./ sum(loopPointWeights,2);
loopPointWeightsNormalized = loopPointWeights ./ max(loopPointWeights,[],2);

[maxWeight, maxLoop] = max(loopPointWeightsNormalized);

bestLoop = zeros(1,size(validWeights,3));
for i = 1:length(trajectoryIndices)
    bestLoop(trajectoryIndices{i}) = finalCommunities(i);
end

bestLoop(bestLoop == 0) = maxLoop(bestLoop == 0);

bestPhases = [];
bestPhase = [];
finalBinPhases = [];
for i = 1:size(loopPhases,1)
%     orderedPhases = sort(loopPhases(i,:));
%     
%     pointsPerBin = ceil(sum(~isnan(orderedPhases)) / NUM_SEGMENTS);
%     nonNanPhases = orderedPhases(~isnan(orderedPhases));
%     
%     binPhases = [nonNanPhases(1:pointsPerBin:end) nonNanPhases(end)];
    
    bestPhases(i,:) = discretize(loopPhases(i,:), loopPhase);
    
%     nanGroups = regionprops(isnan(bestPhases(i,:)), 'PixelIdxList');
%     for j = 1:length(nanGroups)
%         if nanGroups(j).PixelIdxList(1) > 1 && nanGroups(j).PixelIdxList(end) < size(bestPhases,2)
%             bestPhases(i,nanGroups(j).PixelIdxList) = round(linspace(bestPhases(i,nanGroups(j).PixelIdxList(1)-1), bestPhases(i,nanGroups(j).PixelIdxList(end)+1), length(nanGroups(j).PixelIdxList)));
%         end
%     end
    
    stateOccupancy = histcounts(bestPhases(i,:), max(bestPhases(i,:)));
end
for i = 1:size(loopPhases,2)
    bestPhase(i) = bestPhases(bestLoop(i),i);
end




maxPhase = max(bestPhase);
maxLoop = max(bestLoop);
bestStates = (bestLoop-1) * maxPhase + bestPhase; 

index = 1;
stateLoopMap = [];
linearWeights = [];
stateMeans = [];
for i = 1:size(validWeights,1)
    for j = 1:size(validWeights,2)
        linearWeights(index,:) = validWeights(i,j,:);
        stateLoopMap(index,:) = [i,j];
        stateMeans(index,:) = loopMeans{i}(j,:);
        
        index = index + 1;
    end
end
linearWeights = linearWeights ./ sqrt(nansum(linearWeights.^2,2));
% linearWeights = linearWeights ./ max(linearWeights,[],2);
linearWeights(isnan(linearWeights)) = 0;

loopSpaceSimilarities = linearWeights * eye(size(asymmetricProbabilities)) * linearWeights';
% loopSpaceSimilarities = loopSpaceSimilarities ./ sum(loopSpaceSimilarities,2);

loopSpaceOtherSimilarities = loopSpaceSimilarities;
for i = 1:size(validWeights,1)
    thisIDs = find(stateLoopMap(:,1) == i);
    
    loopSpaceOtherSimilarities(thisIDs, thisIDs) = 0;
end

mergeStates = {};
for i = 1:size(loopSpaceOtherSimilarities,1)
    [peaks, indicies] = findpeaks(loopSpaceOtherSimilarities(i,:));
    
    mergeIndicies = indicies(peaks > SIMILARITY_CUTOFF);
    
    [peaks, indicies] = findpeaks(loopSpaceOtherSimilarities(:,i));
    
    mergeIndicies = unique([mergeIndicies, indicies(peaks > SIMILARITY_CUTOFF)]);
    
    mergeStates{i} = mergeIndicies;
end

% oneOffDiagonal = eye(size(asymmetricProbabilities));
% oneOffDiagonal = circshift(oneOffDiagonal,1,2);
% oneOffDiagonal(end,1) = 0;
% 
% loopSpaceTransitions = linearWeights.^2 * oneOffDiagonal * linearWeights.^2';
% loopSpaceTransitions = loopSpaceTransitions ./ sum(loopSpaceTransitions,2);


phaseStates = bestStates;

% [~, betFitStates] = max(linearWeights,[],1);
% bestStates(isnan(bestStates)) = betFitStates(isnan(bestStates));
% bestStates = betFitStates;

% bestLoop = stateLoopMap(bestStates, 1);
% bestPhase = stateLoopMap(bestStates, 2);

finalStates = bestStates;
% for i = 1:length(mergeStates)
%     for j = 1:length(mergeStates{i})
%         finalStates(finalStates == mergeStates{i}(j)) = i;
%     end
% end

stateOccupancy = histcounts(finalStates, max(finalStates));

loopSpaceTransitions = zeros(size(stateLoopMap,1), size(stateLoopMap,1));
for i = 1:size(stateLoopMap,1)    
    if isnan(finalStates(i))
        loopSpaceTransitions(i,i+1) = 1;
        continue;
    end
    
    thisCommunity = stateLoopMap(i,1);
    thisPhase = stateLoopMap(i,2);
    
%     thisIDs = find(linearWeights(i,:) > 0);
    thisIDs = find(finalStates == i);
    
    nextIDs = thisIDs+1;    
    badID = find(nextIDs > length(bestLoop));
    nextIDs(badID) = [];
    thisIDs(badID) = [];
    
    badID = find(~stateHasNext(thisIDs));
    nextIDs(badID) = [];
    thisIDs(badID) = [];    
%     badID = find(~stateHasNext(nextIDs));
%     nextIDs(badID) = [];
%     thisIDs(badID) = [];
    
    badID = find(isnan(finalStates(nextIDs)));
    nextIDs(badID) = [];
    thisIDs(badID) = [];
    
%     badID = find(finalDynamicsStream(thisIDs,end) > 0);
%     nextIDs(badID) = [];
%     thisIDs(badID) = [];
%     badID = find(finalDynamicsStream(nextIDs,end) > 0);
%     nextIDs(badID) = [];
%     thisIDs(badID) = [];
    
    thisWeights = linearWeights(i,thisIDs);
    thisWeights = thisWeights / sum(thisWeights);
    
    nextStates = finalStates(nextIDs);
    
    for j = 1:length(nextStates)
%         loopSpaceTransitions(i,nextStates(j)) = loopSpaceTransitions(i,nextStates(j)) + neighborhoodQuality(thisIDs(j));
%         loopSpaceTransitions(i,nextStates(j)) = loopSpaceTransitions(i,nextStates(j)) + loopPointWeights(thisCommunity,thisIDs(j));
        loopSpaceTransitions(i,nextStates(j)) = loopSpaceTransitions(i,nextStates(j)) + 1;
    end
    
    if nansum(loopSpaceTransitions(i,:)) == 0
        newPhase = mod(i+maxPhase-2, maxPhase)+1;
        newCommunity = thisCommunity;
        loopSpaceTransitions(i,newPhase + (newCommunity-1) * maxPhase) = 1;
    else
        loopSpaceTransitions(i,:) = loopSpaceTransitions(i,:) / nansum(loopSpaceTransitions(i,:));
    end
    if loopSpaceTransitions(i,i) == 1
        loopSpaceTransitions(i,i+1) = 1;
    end
end
loopSpaceTransitions(isnan(loopSpaceTransitions)) = 0;
loopSpaceTransitions = loopSpaceTransitions ./ nansum(loopSpaceTransitions,2);

% figure(2);
% clf;
% hold on;
% h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
% h.Color(4) = 0.2;
% % scatter3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 32, linearWeights(33,:));
% % for i = goodLoops
% %     plot3(loopMeans{i}*pcaBasis(:,1), loopMeans{i}*pcaBasis(:,2), loopMeans{i}*pcaBasis(:,3), 'LineWidth', 2);
% % end
% offsetDistances = [];
% for i = 1:length(bestCommunity)-1
%     if stateHasNext(i)
% %         thisMean = loopMeans{goodLoops(bestCommunity(i))}(bestPhase(i),:);
% %         nextMean = loopMeans{goodLoops(bestCommunity(i+1))}(bestPhase(i+1),:);
% % 
% %         points = [thisMean; nextMean];
%         
%         thisMean = loopMeans{goodLoops(bestCommunity(i))}(bestPhase(i),:);
%         thisPoint = finalDynamicsStream(i,:);
%         
%         points = [thisMean; thisPoint];
%             
%         offsetDistances(i) = sqrt(sum(diff(points*pcaBasis(:,1:3)).^2));
%         
%         plot3(points*pcaBasis(:,1), points*pcaBasis(:,2), points*pcaBasis(:,3), 'k', 'LineWidth', 2);
%     end
% end

simulatedData = finalStates(1);
for t = 2:20000
    simulatedData(t) = randsample(1:size(loopSpaceTransitions,1), 1, true, loopSpaceTransitions(simulatedData(t-1),:));
end

finalSimulatedTrace = stateMeans(simulatedData,:);
discretizedTrace = nan(size(finalDynamicsStream));
discretizedTrace(~isnan(finalStates),:) = stateMeans(finalStates(~isnan(finalStates)),:);
% discretizedTrace(maxWeight < 0.2,:) = nan;
discretizedTrace = discretizedTrace + normrnd(0, 0.2, size(discretizedTrace));
% discretizedTrace(isnan(phaseStates),:) = nan;

figure(21);
clf
subplot(2,1,1);
plot(finalStates);
subplot(2,1,2);
plot(simulatedData(1:length(finalStates)));

plotIndices = 1:length(finalStates);
plotIndices = 1:size(finalSimulatedTrace,1);

figure(22);
clf;
hold on;
h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
h.Color(4) = 0.2;
% plot3(finalSimulatedTrace(plotIndices,:)*pcaBasis(:,1), finalSimulatedTrace(plotIndices,:)*pcaBasis(:,2), finalSimulatedTrace(plotIndices,:)*pcaBasis(:,3));
plot3(discretizedTrace*pcaBasis(:,1), discretizedTrace*pcaBasis(:,2), discretizedTrace*pcaBasis(:,3));
% for i = 1:length(finalStates)-1
%     lineWidth = 0.01;        
%     if ~isnan(finalStates(i)) && ~isnan(finalStates(i+1))
%         lineWidth = sqrt(loopSpaceTransitions(finalStates(i), finalStates(i+1))) * 5;
%     end
%     if lineWidth > 0.01
%         plot3(discretizedTrace(i:i+1,:)*pcaBasis(:,1), discretizedTrace(i:i+1,:)*pcaBasis(:,2), discretizedTrace(i:i+1,:)*pcaBasis(:,3), 'k', 'LineWidth', lineWidth);
%     end
% end
scatter3(stateMeans(104,:)*pcaBasis(:,1), stateMeans(104,:)*pcaBasis(:,2), stateMeans(104,:)*pcaBasis(:,3), 1000, 'kx');
% scatter3(stateMeans(52,:)*pcaBasis(:,1), stateMeans(52,:)*pcaBasis(:,2), stateMeans(52,:)*pcaBasis(:,3), 1000, 'kx');

figure(23);
clf;
hold on;
h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
h.Color(4) = 0.2;
scatter3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 32, bestLoop);

figure(24);
clf;
hold on;
h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
h.Color(4) = 0.2;
scatter3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 32, finalStates);

figure(25);
clf;
hold on;
h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
h.Color(4) = 0.2;
scatter3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 32, isnan(phaseStates));

figure(26);
clf;
hold on;
h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
h.Color(4) = 0.2;
scatter3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 32, loopPhases(3,:));

figure(27);
clf;
hold on;
h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
h.Color(4) = 0.2;
plot3(finalSimulatedTrace(plotIndices,:)*pcaBasis(:,1), finalSimulatedTrace(plotIndices,:)*pcaBasis(:,2), finalSimulatedTrace(plotIndices,:)*pcaBasis(:,3));

figure(28);
clf;
subplot(2,1,1);
imagesc(finalDynamicsStream');
originalCAxis = caxis;
subplot(2,1,2);
imagesc(finalSimulatedTrace(1:size(finalDynamicsStream,1),:)');
caxis(originalCAxis);

%% Assign inputs

inputChance = [];
for i = 1:size(stateLoopMap,1)
    thisIDs = find(finalStates == i);
    
    inputChance(i) = mean(finalDynamicsStream(thisIDs,end));
end
inputChance(isnan(inputChance)) = 0;

figure(1);
clf;
hold on;
plot(inputChance);
plot(stateOccupancy/max(stateOccupancy))

figure(2);
clf;
hold on;
h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
h.Color(4) = 0.2;
for i = 1:length(loopMeans)
    thisIDs = find(stateLoopMap(:,1) == i);
    
%     sizes = stateOccupancy(thisIDs)/max(stateOccupancy(thisIDs));
%     sizes = inputChance(thisIDs);
    sizes = 1:length(thisIDs);
    sizes = sizes / max(sizes) * 100;
    sizes = sizes + 32;
%     scatter3(loopMeans{i}(1:end-1,:)*pcaBasis(:,1), loopMeans{i}(1:end-1,:)*pcaBasis(:,2), loopMeans{i}(1:end-1,:)*pcaBasis(:,3), sizes, inputChance(thisIDs));
    scatter3(loopMeans{i}(1:end-1,:)*pcaBasis(:,1), loopMeans{i}(1:end-1,:)*pcaBasis(:,2), loopMeans{i}(1:end-1,:)*pcaBasis(:,3), sizes, sizes);
end

%%

% figure(5);
% clf;
% hold on;
% h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
% h.Color(4) = 0.2;
% colorValues = loopCounts;
% % colorValues(colorValues < -10) = -10;
% caxis([min(colorValues) max(colorValues)]);
% colorValues = colorValues - min(colorValues);
% colorValues = colorValues / max(colorValues);
% colorValues = round(colorValues*255 + 1);
% colors = jet(256);
% colormap(jet(256));
% colorbar
% for i = 37%1:length(trajectoryMeans)
% %     if sum(mergedTrajectories(i,:)) + sum(mergedTrajectories(:,i)) <= 0
% %         continue;
% %     end
%
%
%     maxValue = max(mergedStartConnections(i,:));
%     startTrajectories = find(mergedStartConnections(i,:) / maxValue > CONNECT_CUTOFF & mergedTrajectories(i,:));
%
%     maxValue = max(mergedEndConnections(i,:));
%     endTrajectories = find(mergedEndConnections(i,:) / maxValue > CONNECT_CUTOFF & mergedTrajectories(i,:));
%
%     for j = 1:length(startTrajectories)
%         points(1,:) = trajectoryMeans(1,:,i);
%         points(2,:) = trajectoryMeans(end,:,startTrajectories(j));
%
%         plot3(points*pcaBasis(:,1), points*pcaBasis(:,2), points*pcaBasis(:,3), 'k', 'LineWidth', 1);
%     end
%
%     for j = 1:length(endTrajectories)
%         points(1,:) = trajectoryMeans(end,:,i);
%         points(2,:) = trajectoryMeans(1,:,endTrajectories(j));
%
%         plot3(points*pcaBasis(:,1), points*pcaBasis(:,2), points*pcaBasis(:,3), 'k', 'LineWidth', 1);
%     end
%
%     plot3(trajectoryMeans(:,:,i)*pcaBasis(:,1), trajectoryMeans(:,:,i)*pcaBasis(:,2), trajectoryMeans(:,:,i)*pcaBasis(:,3), 'Color', colors(colorValues(i),:), 'LineWidth', 3);
%
% end
%
% %%
%
% figure(5);
% clf;
% hold on;
% h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
% h.Color(4) = 0.2;
% colorValues = (connectionTrajectoryCount);
% colorValues = colorValues - min(colorValues(:));
% colorValues = colorValues / max(colorValues(:));
% colorValues = round(colorValues*255 + 1);
% colors = jet(256);
% colormap(jet(256));
% caxis([1 256]);
% colorbar
%
% for i = 1:length(trajectoryMeans)
%     thisIDs = find(connectionTrajectoryCount(i,:) > 1);
%
%     for j = 1:length(thisIDs)
%         thisTrajectory = thisIDs(j);
%
%         if thisTrajectory == i
%             continue;
%         end
%
%         if endTrajectoryCount(j,i) > startTrajectoryCount(i,j)
%             points(1,:) = trajectoryMeans(1,:,i);
%             points(2,:) = trajectoryMeans(end,:,thisTrajectory);
%         else
%             points(1,:) = trajectoryMeans(1,:,i);
%             points(2,:) = trajectoryMeans(end,:,thisTrajectory);
%         end
%
%         color = colors(colorValues(i,thisTrajectory),:);
%
%         plot3(points*pcaBasis(:,1), points*pcaBasis(:,2), points*pcaBasis(:,3), 'Color', color, 'LineWidth', 3);
%     end
%
%     plot3(trajectoryMeans(:,:,i)*pcaBasis(:,1), trajectoryMeans(:,:,i)*pcaBasis(:,2), trajectoryMeans(:,:,i)*pcaBasis(:,3), 'k', 'LineWidth', 2);
% end
%
% %%
%
% CONNECT_CUTOFF = 0.01;
%
% startConnectedTrajectories = [];
% endConnectedTrajectories = [];
% for i = 1:size(mergedStartConnections,1)
%     maxValue = max(mergedStartConnections(i,:));
%     startConnectedTrajectories(i,:) = mergedStartConnections(i,:) / maxValue > CONNECT_CUTOFF;
%
%     maxValue = max(mergedEndConnections(i,:));
%     endConnectedTrajectories(i,:) = mergedEndConnections(i,:) / maxValue > CONNECT_CUTOFF;
% end
% mergedTrajectories = min(startConnectedTrajectories, endConnectedTrajectories');
%
% figure(2);
% clf;
% colormap(jet(256));
% subplot(2,2,1);
% imagesc(startConnectedTrajectories);
% title("Start connections");
% subplot(2,2,2);
% imagesc(endConnectedTrajectories');
% title("End connections");
% subplot(2,2,3);
% imagesc(mergedTrajectories');
% title("Merge connections");
%
% connectionQualities = [];
% totalQuality = [];
% for i = 1:length(trajectoryMeans)
% %     connectionQualities(i) = sum(mergedEndConnections(i,:) .* trajectoryCounts) .* sum(mergedStartConnections(i,:) .* trajectoryCounts);
% %     totalQuality(i) = connectionQualities(i) .* trajectoryCounts(i);
%     maxStartValue = max(mergedStartConnections(i,:));
%     maxEndValue = max(mergedEndConnections(i,:));
%     connectionQualities(i) = mean(mergedEndConnections(i,:) / maxEndValue .* trajectoryCounts) + mean(mergedStartConnections(i,:) / maxStartValue .* trajectoryCounts);
%
%     totalQuality(i) = connectionQualities(i) .* trajectoryCounts(i);
% %     totalQuality(i) = sqrt(connectionQualities(i).^2 + (2*trajectoryCounts(i)).^2);
% end
%
% figure(5);
% clf;
% hold on;
% h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
% h.Color(4) = 0.2;
% colorValues = log(totalQuality);
% colorValues = colorValues - min(colorValues);
% colorValues = colorValues / max(colorValues);
% colorValues = round(colorValues*255 + 1);
% colors = jet(256);
% colormap(jet(256));
% caxis([1 256]);
% colorbar
% for i = 1:length(trajectoryMeans)
% %     if sum(mergedTrajectories(i,:)) + sum(mergedTrajectories(:,i)) <= 0
% %         continue;
% %     end
%
%
%     maxValue = max(mergedStartConnections(i,:));
%     startTrajectories = find(mergedStartConnections(i,:) / maxValue > CONNECT_CUTOFF & mergedTrajectories(i,:));
%
%     maxValue = max(mergedEndConnections(i,:));
%     endTrajectories = find(mergedEndConnections(i,:) / maxValue > CONNECT_CUTOFF & mergedTrajectories(i,:));
%
%     for j = 1:length(startTrajectories)
%         points(1,:) = trajectoryMeans(1,:,i);
%         points(2,:) = trajectoryMeans(end,:,startTrajectories(j));
%
%         plot3(points*pcaBasis(:,1), points*pcaBasis(:,2), points*pcaBasis(:,3), 'k', 'LineWidth', 1);
%     end
%
%     for j = 1:length(endTrajectories)
%         points(1,:) = trajectoryMeans(end,:,i);
%         points(2,:) = trajectoryMeans(1,:,endTrajectories(j));
%
%         plot3(points*pcaBasis(:,1), points*pcaBasis(:,2), points*pcaBasis(:,3), 'k', 'LineWidth', 1);
%     end
%
%     plot3(trajectoryMeans(:,:,i)*pcaBasis(:,1), trajectoryMeans(:,:,i)*pcaBasis(:,2), trajectoryMeans(:,:,i)*pcaBasis(:,3), 'Color', colors(colorValues(i),:), 'LineWidth', 3);
%
% end
%
% %%
%
% [steadyState, ~] = eigs(asymmetricProbabilities', 1);
% fluxMatrix = diag(steadyState) * asymmetricProbabilities;
% symmetricMatrix = (fluxMatrix+fluxMatrix.')/2;
% antisymmetricMatrix = (fluxMatrix-fluxMatrix.')/2;
% symmetricMatrix = symmetricMatrix ./ steadyState;
% antisymmetricMatrix = antisymmetricMatrix ./ steadyState;
%
% [eigenVectors,~] = eigs(symmetricMatrix,10);
% plotDistribution = eigenVectors(:,2);
%
% figure(6);
% clf;
% hold on;
% h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
% h.Color(4) = 0.2;
% scatter3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 32, plotDistribution);
%
%
% %%
%
%
% % weightOverlap = stateWeights * stateWeights';
% % weightOverlap = weightOverlap^2;
%
% % KLWeights = stateWeights;
% % KLWeights(KLWeights < eps) = eps;
% % waitHandle = parfor_progressbar(size(KLWeights, 1), 'Calculating KLs');
% % weightOverlap = [];
% % for i = 1:size(KLWeights, 1)
% %     klValues = real(log(KLWeights ./ KLWeights(i,:)));
% %     klValues(isinf(klValues)) = 0;
% %     klValues(isnan(klValues)) = 0;
% %
% %     weightOverlap(i,:) = -sum(KLWeights(i,:) .* klValues,2);
% %
% %     waitHandle.iterate(1);
% % end
% % close(waitHandle);
% %
% % weightOverlap = weightOverlap .* weightOverlap';
%
% %%
% %
% % figure(4);
% % clf;
% % hold on;
% % h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
% % h.Color(4) = 0.2;
% % colors = jet(length(trajectoryMeans));
% % for i = 35%1:length(trajectoryMeans)
% %     meanValues = trajectoryMeans(:,:,i);
% %     plot3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 'LineWidth', 3, 'Color', colors(i,:));
% % %     scatter3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 32, trajectoryQuality(:,i));
% % %         scatter3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 32, trajectoryCertainty(:,i));
% %
% %     thisStateIDs = trajectoryStateMap(i,:);
% %
% %     for j = 1:size(trajectoryMeans,1)
% %         stateID = trajectoryStateMap(i,j);
% %
% %         transitionStates = weightOverlap(stateID,:);
% %         maxValue = max(transitionStates);
% %         transitionStates(thisStateIDs) = 0;
% %         transitionStates = find(transitionStates / maxValue > 0.75);
% %
% %         transitionStates = transitionStates(1:10:end);
% %
% %         for k = 1:length(transitionStates)
% %             tempIDs = stateTrajectoryMap(transitionStates(k),:);
% %             toMeans = squeeze(trajectoryMeans(tempIDs(2), :, tempIDs(1)));
% %             points = [squeeze(meanValues(j,:)); toMeans];
% %             plot3(points*pcaBasis(:,1), points*pcaBasis(:,2), points*pcaBasis(:,3), 'k', 'LineWidth', 0.5);
% %         end
% %     end
% % end
% % colormap(jet(256));
% % caxis([1 length(trajectoryMeans)])
% %
% %% Link clusters
%
% normWeights = stateWeights ./ sqrt(sum(stateWeights.^2,2));
% weightOverlap = normWeights * normWeights';
% % weightOverlap = weightOverlap^2;
%
% % symmetricSimilarities = min(stateSimilarities,stateSimilarities');
% % outOfClusterTransitions = [];
% outOfClusterSimilarities = zeros(size(stateSimilarities));
% % stateOverlap = [];
%
% figure(5);
% clf;
% hold on;
% h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
% h.Color(4) = 0.2;
% for i = 1:length(trajectoryMeans)
%     meanValues = trajectoryMeans(:,:,i);
%
%     thisStateIDs = trajectoryStateMap(i,:);
%
%     for j = 1:size(trajectoryMeans,1)
%         stateID = trajectoryStateMap(i,j);
%
%         transitionStates = weightOverlap(stateID,:);
%         transitionStates = [0 transitionStates 0];
% %         maxValue = max(transitionStates);
% %         transitionStates(thisStateIDs) = 0;
% %         outOfClusterSimilarities(stateID, :) = transitionStates;
%         [peakValues, transitionStatesIDs] = findpeaks(transitionStates);
%         transitionStatesIDs = transitionStatesIDs - 1;
%         transitionStatesIDs = transitionStatesIDs(find(peakValues > 0.1));
%         transitionStates = transitionStates(transitionStatesIDs);
%         outOfClusterSimilarities(stateID, transitionStatesIDs) = 1;
% %         transitionStatesIDs = 1:length(transitionStates);
% %         transitionStates(thisStateIDs) = transitionStates(thisStateIDs) .* trajectoryQuality(j,i)^2;
%         transitionStatesWeights = transitionStates ./ sum(transitionStates);
%
% %         outOfClusterTransitions(stateID) = sum(transitionStates);
%
%         transitionStateMeans = zeros(size(meanValues(j,:)));
%         for k = 1:length(transitionStates)
%             transitionStateID = stateTrajectoryMap(transitionStatesIDs(k),:);
%             thisMean = squeeze(trajectoryMeans(transitionStateID(2),:,transitionStateID(1)));
%
%             transitionStateMeans = transitionStateMeans + transitionStatesWeights(k).*thisMean;
%         end
%
% %         meanValues(j,:) = meanValues(j,:) * trajectoryQuality(j,i) + (1 - trajectoryQuality(j,i)) * transitionStateMeans;
%         meanValues(j,:) = transitionStateMeans;
%     end
%
%      plot3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 'LineWidth', 3);
%      plot3(trajectoryMeans(:,:,i)*pcaBasis(:,1), trajectoryMeans(:,:,i)*pcaBasis(:,2), trajectoryMeans(:,:,i)*pcaBasis(:,3), 'k', 'LineWidth', 1);
%
% end
% %
% % figure(6);
% % clf
% % imagesc(outOfClusterSimilarities);
% %


%%

% clusterMask = kron(eye(length(trajectoryMeans)), ones(size(trajectoryMeans,1)));
% stateDistances = ~clusterMask .* stateSimilarities;
% % stateDistances(stateDistances == 0) = max(stateDistances(:));
% stateDistances(stateDistances < eps) = eps;
% stateDistances = 1./stateDistances;
% stateDistances = stateDistances .* ~eye(size(stateDistances));
% stateDistances = squareform(stateDistances);
%
% clusteringData = linkage(stateDistances, 'complete');
%
% figure(1)
% clf
% plot(((clusteringData(1:500,3))))
%
% clusterIDs = cluster(clusteringData, 'cutoff', 200,'Criterion','distance');
%
% %%
%
% reducedMeans = [];
% for i = 1:max(clusterIDs)
%     thisStates = find(clusterIDs == i);
%
%     thisWeights = 1./weightOverlap(thisStates, thisStates);
%     thisWeights(isinf(thisWeights)) = 1;
%     thisWeights = sum(thisWeights);
%     thisWeights = thisWeights .* (trajectoryCounts(stateTrajectoryMap(thisStates,2))).^2;
% %     thisWeights = ones(size(thisWeights));
%     thisWeights = thisWeights ./ sum(thisWeights);
%
%     meanWeights = thisWeights * stateWeights(thisStates, :);
%     meanWeights = meanWeights / sum(meanWeights);
% %     meanWeights = stateWeights(thisStates(1), :);
%
% %     thisReducedMean = zeros(size(trajectoryMeans(1,:,1)));
% %     for k = 1:length(thisStates)
% %         transitionStateID = stateTrajectoryMap(thisStates(k),:);
% %         thisMean = squeeze(trajectoryMeans(transitionStateID(2),:,transitionStateID(1)));
% %
% %         thisReducedMean = thisReducedMean + thisWeights(k).*thisMean;
% %     end
%
%     reducedMeans(i,:) = meanWeights * finalDynamicsStream;
%
%     stateTrajectories = stateTrajectoryMap(thisStates,:);
%     stateMeans = [];
%     for j = 1:size(stateTrajectories,1)
%         stateMeans(j,:) = trajectoryMeans(stateTrajectories(j,2),:,stateTrajectories(j,1));
%     end
%
%     if length(thisStates) > 1 && sum(stateTrajectories(:,1) == 35) > 0
%         figure(5);
%         clf;
%         hold on;
%         h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
%         h.Color(4) = 0.2;
%         scatter3(stateMeans*pcaBasis(:,1), stateMeans*pcaBasis(:,2), stateMeans*pcaBasis(:,3), 'x');
%         scatter3(reducedMeans(end,:)*pcaBasis(:,1), reducedMeans(end,:)*pcaBasis(:,2), reducedMeans(end,:)*pcaBasis(:,3));
%     end
% end
%
% figure(5);
% clf;
% hold on;
% h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
% h.Color(4) = 0.2;
% scatter3(reducedMeans*pcaBasis(:,1), reducedMeans*pcaBasis(:,2), reducedMeans*pcaBasis(:,3));

%%


%% Merge points


%
% clusterSimilarities = outOfClusterSimilarities;
% clusterSimilarities = min(clusterSimilarities, clusterSimilarities');
% weightSimilarities = weightOverlap;
%
% oldStateMap = {};
% stateToNewMap = [];
%
% for i = 1:size(weightOverlap,1)
%     oldStateMap{i} = i;
% end
%
%
% while 1
%     bestStateLists = {};
%     statesLeft = 1:size(weightSimilarities,1);
%     waitHandle = parfor_progressbar(size(weightSimilarities,1), 'Merging trajectories');
%     while length(statesLeft) > 0
%         nextStateID = statesLeft(1);
%
%         stateLists = {};
%         stateListValues = [];
%
%         thisID = nextStateID;
%         currentStateList = thisID;
%         checkedList = thisID;
%
%         outStates = find(clusterSimilarities(currentStateList,:) == 1);
%         inStates = find(clusterSimilarities(:,currentStateList) == 1);
%         outStates = intersect(outStates, statesLeft);
%         inStates = intersect(inStates, statesLeft);
%
%         nextStateList = unique([currentStateList, outStates,  inStates']);
%         diffStateList = setdiff(nextStateList, currentStateList);
%
%         stateLists{1} = currentStateList;
%         inClusterValues = weightSimilarities(currentStateList, currentStateList);
%         outClusterValues = weightSimilarities(currentStateList, diffStateList);
%         stateListValues(1) = sum(inClusterValues(:)) / (sum(inClusterValues(:)) + sum(outClusterValues(:)));
%
%         uncheckedStateList = setdiff(nextStateList, checkedList);
%
%         while length(uncheckedStateList) > 0 && length(currentStateList) <= 10
%             clusterWeights = weightSimilarities(currentStateList, uncheckedStateList);
%             clusterWeights = sum(clusterWeights);
%             [~, thisID] = max(clusterWeights);
%             thisID = uncheckedStateList(thisID);
%             currentStateList = [currentStateList thisID];
%             checkedList = [checkedList thisID];
%
%             outStates = [];
%             inStates = [];
%             for j = 1:length(currentStateList)
%                 outStates = [outStates find(clusterSimilarities(currentStateList(j),:) == 1)];
%                 inStates = [inStates; find(clusterSimilarities(:,currentStateList(j)) == 1)];
%             end
%             outStates = intersect(unique(outStates), statesLeft);
%             inStates = intersect(unique(inStates), statesLeft);
%
%             nextStateList = unique([currentStateList, outStates,  inStates']);
%             diffStateList = setdiff(nextStateList, currentStateList);
%
%             stateLists{end+1} = currentStateList;
%             inClusterValues = weightSimilarities(currentStateList, currentStateList);
%             outClusterValues = weightSimilarities(currentStateList, diffStateList);
%             stateListValues(end+1) = sum(inClusterValues(:)) / (sum(inClusterValues(:)) + sum(outClusterValues(:)));
%
%             uncheckedStateList = setdiff(nextStateList, checkedList);
%         end
%
%         [~, bestListID] = max(stateListValues);
%         bestStateList = [];
%         for i = 1:length(stateLists{bestListID})
%             bestStateList = [bestStateList oldStateMap{stateLists{bestListID}(i)}];
%         end
%         bestStateLists{end+1} = bestStateList;
%         stateToNewMap(bestStateLists{end}) = length(bestStateLists);
%
%         statesLeft = setdiff(statesLeft, stateLists{bestListID});
%         waitHandle.iterate(length(stateLists{bestListID}));
%     end
%     close(waitHandle);
%
%     newClusterSimilarities = zeros(length(bestStateLists), length(bestStateLists));
%     newWeightSimilarities = weightOverlap;
%     removedIndices = [];
%     for i = 1:length(bestStateLists)
%         mergeSimilarities = find(sum(outOfClusterSimilarities(bestStateLists{i},:), 1) >= 1);
%         mergeSimilarities = unique(stateToNewMap(mergeSimilarities));
%
%         firstIndex = bestStateLists{i}(1);
%         newWeightSimilarities(firstIndex,:) = mean(newWeightSimilarities(bestStateLists{i},:), 1);
%         newWeightSimilarities(:,firstIndex) = mean(newWeightSimilarities(:,bestStateLists{i}), 2);
%         if length(bestStateLists{i}) > 1
%             removedIndices = [removedIndices bestStateLists{i}(2:end)];
%         end
%
%         newClusterSimilarities(i, mergeSimilarities) = 1;
%     end
%     newWeightSimilarities(removedIndices,:) = [];
%     newWeightSimilarities(:,removedIndices) = [];
%
%     if size(newClusterSimilarities,1) == size(clusterSimilarities,1)
%         break;
%     end
%
%     clusterSimilarities = newClusterSimilarities;
%     weightSimilarities = newWeightSimilarities;
%     oldStateMap = bestStateLists;
%
%     disp(["Clusters : " + size(clusterSimilarities,1) + " weights: " + size(weightSimilarities,1) + " next step: " + length(oldStateMap)]);
% end
%
% %%
%
% reducedMeans = [];
% for i = 1:length(bestStateLists)
%     thisWeights = weightOverlap(bestStateLists{i}, bestStateLists{i});
%     thisWeights = sum(thisWeights);
%     thisWeights = ones(size(thisWeights));
%     thisWeights = thisWeights ./ sum(thisWeights);
%
%     thisReducedMean = zeros(size(trajectoryMeans(1,:,1)));
%     for k = 1:length(bestStateLists{i})
%         transitionStateID = stateTrajectoryMap(bestStateLists{i}(k),:);
%         thisMean = squeeze(trajectoryMeans(transitionStateID(2),:,transitionStateID(1)));
%
%         thisReducedMean = thisReducedMean + thisWeights(k).*thisMean;
%     end
%
%     reducedMeans(i,:) = thisReducedMean;
% end
%
% figure(5);
% clf;
% hold on;
% h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
% h.Color(4) = 0.2;
% scatter3(reducedMeans*pcaBasis(:,1), reducedMeans*pcaBasis(:,2), reducedMeans*pcaBasis(:,3));
%
% figure(6)
% clf;
% hold on
% h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
% h.Color(4) = 0.2;
% % colormap(lines(256));
% % colormap(jet(256));
% colormap(hsv(256));
% for i = 1:length(bestStateLists)
%     thisPoints = [];
%     phases = [];
%     for k = 1:length(bestStateLists{i})
%         transitionStateID = stateTrajectoryMap(bestStateLists{i}(k),:);
%         thisPoints(k,:) = squeeze(trajectoryMeans(transitionStateID(2),:,transitionStateID(1)));
%
%         thisWeights = squeeze(trajectoryWeights(:,transitionStateID(2),transitionStateID(1)));
%         phases(k) = angle(sum(thisWeights' .* exp(1i*allPhases)));
%     end
%
%     scatter3(thisPoints*pcaBasis(:,1), thisPoints*pcaBasis(:,2), thisPoints*pcaBasis(:,3), 32, phases);
% %     scatter3(thisPoints*pcaBasis(:,1), thisPoints*pcaBasis(:,2), thisPoints*pcaBasis(:,3), 32, ones(size(phases))'*i);
% end
%
% % rawCommunities = community_louvain(outOfClusterSimilarities, 1000);
% % [~, sortIDs] = sort(rawCommunities);
% % imagesc(outOfClusterSimilarities(sortIDs,sortIDs));
%
% %%
%
% stateOverlap = max(stateOverlap, stateOverlap');
% stateOverlap(stateOverlap < 0) = 0;
%
% figure(6);
% clf
% % hold on;
% % plot(outOfClusterTransitions)
% % plot(outOfClusterTransitions)
% imagesc(weightOverlap);
%
% rawCommunities = community_louvain(weightOverlap^5, 10);
% [~, sortIDs] = sort(rawCommunities);
% imagesc(weightOverlap(sortIDs,sortIDs));
%
% mergeStates = find(max(outOfClusterSimilarities, [], 2)' .* (1 - outOfClusterTransitions) > 0.1);
% % mergeStates = find(outOfClusterTransitions > 0.5);
%
% figure(7);
% clf;
% hold on;
% h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
% h.Color(4) = 0.2;
% colors = lines(max(rawCommunities));
% for i = 1:length(trajectoryMeans)
%     meanValues = trajectoryMeans(:,:,i);
%
%     thisStateIDs = trajectoryStateMap(i,:);
%     thisStates = 1:length(thisStateIDs);%find(ismember(thisStateIDs, mergeStates));
%
%     plot3(meanValues*pcaBasis(:,1), meanValues*pcaBasis(:,2), meanValues*pcaBasis(:,3), 'k', 'LineWidth', 3);
%     scatter3(meanValues(thisStates,:)*pcaBasis(:,1), meanValues(thisStates,:)*pcaBasis(:,2), meanValues(thisStates,:)*pcaBasis(:,3), 32, colors(rawCommunities(thisStateIDs(thisStates)),:));
%
% end
%
% %%
% mergeStates = zeros(size(stateOverlap));
%
% for i = 1:length(trajectoryMeans)
%     thisStateIDs = trajectoryStateMap(i,:);
%     for j = i+1:length(trajectoryMeans)
%         transitionStateIDs = trajectoryStateMap(j,:);
%
%         thisTransitions = stateOverlap(thisStateIDs,transitionStateIDs);
%
% %         while max(thisTransitions(:)) > 0.01
% %             [~, maxIndex] = max(thisTransitions,[],2);
% %             [maxX, maxY] = ind2sub(size(thisTransitions), maxIndex);
% %
% %             thisTransitions(maxX, maxY);
% %         end
%
%         if max(thisTransitions(:)) > 0.01
%             figure(1);
%             clf;
%             imagesc(thisTransitions)
%             caxis([min(stateOverlap(:)) max(stateOverlap(:))])
%             colorbar
%         end
%     end
% end

