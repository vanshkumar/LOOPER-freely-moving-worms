
 
if ~exist('clusterCounts', 'var') || isempty(clusterCounts)
    clusterCounts = [100 80 60 50 40 30 20 17 15 12 10 8];
end

if ~exist('distanceType', 'var') || isempty(distanceType)
    distanceType = 'correlation';
end

if ~exist('maxCheckTime', 'var') || isempty(maxCheckTime)
    maxCheckTime = 10;
end

hadError = 0;

try
    pdist2(1,1,distanceType);
catch
    msgbox('Distance measure must be one of: euclidean squaredeuclidean seuclidean mahalanobis cityblock minkowski chebychev cosine correlation hamming jaccard or spearman. See pdist2 documentation for details.');
    
    hadError = 1;
    
    return;
end


%%

showIndices = 1:size(finalDynamicsStream,1);
FIG_BASE = 100;


diffusedProbabilities = asymmetricProbabilities;

%% Get time evolution of true matrix

MIN_MODE = 0.1;
MAX_TIMES = 50;
PHASE_SIGMA = pi/10;
PHASE_STEPS = 100;
checkTimes = 1:maxCheckTime;

stepMatrix = diffusedProbabilities;

stateDistributions = [];
currentMatrix = diffusedProbabilities^checkTimes(1);

waitHandle = parfor_progressbar(length(checkTimes), 'Calculating true distributions');
for i = 1:length(checkTimes)
    stateDistributions(:,:,i) = currentMatrix*eye(size(diffusedProbabilities,2));
    
    currentMatrix = stepMatrix*currentMatrix;
    
    waitHandle.iterate(1);
end
close(waitHandle);




%% Build similarity of state matrix

allTimeMatrix = (asymmetricProbabilities);
similarities1 = pdist((allTimeMatrix), distanceType);

stateDistances = zeros(1, size(diffusedProbabilities,1) * (size(diffusedProbabilities,1)-1) / 2);
statePairIndex = 1;
waitHandle = parfor_progressbar(size(diffusedProbabilities,1), 'Calculating distances');
for i = 1:size(diffusedProbabilities,1)
    
    for j = i+1:size(diffusedProbabilities,1)
    
        stateDistances(statePairIndex) = similarities1(statePairIndex); 
        statePairIndex = statePairIndex + 1;
    end
    waitHandle.iterate(1);
end
close(waitHandle);

clustering = linkage(stateDistances, 'average');


KLDivergences = [];
reducedMatrices = {};
cosineDistances = [];
waitHandle = parfor_progressbar(length(clusterCounts), 'Calculating reduced distributions');
for clusterIndex = 1:length(clusterCounts)
    
    maxClusters = clusterCounts(clusterIndex);
    clusterIDs = cluster(clustering,'maxclust',maxClusters);
    
    
    disp(['Calculating response of ' num2str(maxClusters) ' clusters']);
    
    %% Build reduced matrix
    
    reducedMatrixTemp = zeros(max(clusterIDs), size(asymmetricProbabilities,2));
    for i = 1:size(reducedMatrixTemp,1)
        reducedMatrixTemp(i,:) = sum(asymmetricProbabilities(clusterIDs == i, :), 1) ./ sum(clusterIDs == i);
    end
    
    reducedMatrix = zeros(max(clusterIDs));
    for i = 1:size(reducedMatrix,1)
        reducedMatrix(:,i) = sum(reducedMatrixTemp(:, clusterIDs == i), 2);
    end
    
    reducedStateCounts = [];
    for j = 1:size(reducedMatrix,2)
        reducedStateCounts(j) = sum(clusterIDs == j);
    end
    
    testStateDistributions = [];
    currentMatrix = reducedMatrix^checkTimes(1);
    stepMatrix = reducedMatrix;
    
    for i = 1:length(checkTimes)
        reducedDistribution = currentMatrix*eye(size(reducedMatrix,2));
        
        expandedMatrixTemp = zeros(length(clusterIDs), size(reducedDistribution,2));
        for j = 1:size(expandedMatrixTemp,2)
            numClusters = sum(clusterIDs == j);
            expandedMatrixTemp(clusterIDs == j,:) = repmat(reducedDistribution(j, :), [numClusters,1]);
        end
        
        expandedMatrix = zeros(length(clusterIDs), length(clusterIDs));
        for j = 1:size(expandedMatrixTemp,2)
            numClusters = sum(clusterIDs == j);
            expandedMatrix(:,clusterIDs == j) = repmat(expandedMatrixTemp(:, j), [1,numClusters]) / numClusters;
        end
        
        testStateDistributions(:,:,i) = expandedMatrix;
        
        
        
        currentMatrix = stepMatrix*currentMatrix;
    end
    
    
    %% Calculate KLs
    
    
    for i = 1:length(checkTimes)
        realDistribution = stateDistributions(:,:,i);
        reducedDistribution = testStateDistributions(:,:,i);
        
        elementwiseKL = reducedDistribution.*log(reducedDistribution ./ realDistribution);
        elementwiseKL(isnan(elementwiseKL)) = 0;
        elementwiseKL(isinf(elementwiseKL)) = 0;
        
        
        KLDivergences(clusterIndex,i) = quantile(sum(elementwiseKL,2), 0.95);
    end
    
    figure(FIG_BASE + 1);
    clf;
    plot(KLDivergences(clusterIndex,:));
    title(['KL divergence (clusters = ' num2str(maxClusters) ')']);
    xlabel('Model time step');
    ylabel('95% quantile KL');
    
    
    
    reducedMatrices{clusterIndex} = reducedMatrix;
    waitHandle.iterate(1);
end
close(waitHandle);

%%

worstDivervenges = [];
KLPeaks = [];
for i = 1:size(KLDivergences,1)
    worstDivervenges(i) = max(KLDivergences(i,:));

end

BICs = worstDivervenges*size(diffusedProbabilities,1) + (clusterCounts).^2/2*log(size(diffusedProbabilities,1)/(2*pi));

figure(FIG_BASE + 1);
clf;
plot(KLDivergences');
strings = {};
for i = 1:length(clusterCounts)
    strings{i} = [num2str(clusterCounts(i)) ' clusters'];
end
legend(strings);
title('Time series of model fits (KL divergence)');
xlabel('Model time steps');
ylabel('95% qunatile of KL divergence');

figure(FIG_BASE + 2);
clf
plot(clusterCounts, BICs);
title('Model fitting (BIC vs cluster count)');
xlabel('Number of clusters');
ylabel('BIC');

%% Build best model
[~, minIndex] = min(BICs);
maxClusters = clusterCounts(minIndex);
clusterIndex = find(clusterCounts == maxClusters);
clusterIDs = cluster(clustering,'maxclust',maxClusters);

uniqueClusterIDs = unique(clusterIDs);
finalClusterCount = length(uniqueClusterIDs);

for i = 1:length(uniqueClusterIDs)
    clusterIDs(clusterIDs == uniqueClusterIDs(i)) = i;
end

validClusters = [];
countClusters = [];
clusterMeans = [];
clusterMeansPCA = [];
clusterSTDs = [];
clusterSTDsPCA = [];

maxDims = min(3, size(pcaBasis,2));

for i = 1:max(clusterIDs)
    thisIndices = find(clusterIDs == i);
    
    clusterMeans(i,:) = mean(finalDynamicsStream(thisIndices,:),1);
    clusterSTDs(i,:) = std(finalDynamicsStream(thisIndices,:), [], 1) / maxDims;
    
    clusterMeansPCA(i,:) = mean(finalDynamicsStream(thisIndices,:)*pcaBasis(:,1:maxDims),1);
    clusterSTDsPCA(i,:) = std(finalDynamicsStream(thisIndices,:)*pcaBasis(:,1:maxDims),[],1) / maxDims;
    
    countClusters(i) = length(thisIndices);
    validClusters(i) = countClusters(i) > size(diffusedProbabilities,1) / size(reducedMatrices{clusterIndex},1) / 10;
    
    
    DEBUG = 0;
    if DEBUG
        figure(FIG_BASE + i);
        clf;
        hold on;
        thisTrace = finalDynamicsStream * pcaBasis;
        plot3(thisTrace(:,1), thisTrace(:,2), thisTrace(:,3));
        scatter3(thisTrace(thisIndices,1), thisTrace(thisIndices,2), thisTrace(thisIndices,3));
        scatter3(clusterMeansPCA(i,1), clusterMeansPCA(i,2), clusterMeansPCA(i,3), 300, 'rx', 'LineWidth', 10);
    end
    
end

finalReducedMatrix = reducedMatrices{clusterIndex}(uniqueClusterIDs,uniqueClusterIDs);
