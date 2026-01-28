
%% get variables

if ~exist('repopulateDensity', 'var') || isempty(repopulateDensity)
    repopulateDensity = 0.95;
end
if ~exist('minReturnTime', 'var') || isempty(minReturnTime)
    minReturnTime = 10;
end

FIG_BASE = 200;

    
    %%

    thisDynamicsStream = timeSeriesData';

    if size(thisDynamicsStream,2) >= 3
        [pcaBasis, ~] = pca(thisDynamicsStream, 'NumComponents', 3);
    else
        pcaBasis = eye(size(thisDynamicsStream,2));
    end

    meanDistanceMatrix = pdist2(thisDynamicsStream, thisDynamicsStream);

    bestSigmaValues = [];
    bestNeighborCounts = [];
    bestLocalProjections = [];
    distanceRatios = [];
    bestDiffSigmaValues = [];
    localDistances = [];


    TRAJECTORY_SIZE = 2;        
    CALC_ALL = 1;        

    if ~CALC_ALL
        calcIndicies = 5;
    else
        calcIndicies = 1:size(thisDynamicsStream,1);
    end

    bestDistances = ones(size(calcIndicies))*1/100;
    bestNeighborhood = zeros(size(calcIndicies));
    waitHandle = parfor_progressbar(size(thisDynamicsStream,1), 'Calculating nearest neighbors');        
    parfor calculateIndex = calcIndicies
        [bestLocalProjections(:,calculateIndex), bestSigmaValues(calculateIndex), bestDiffSigmaValues(calculateIndex), localDistances(:,calculateIndex), bestNeighborCounts(calculateIndex), distances] = findBestAxes(thisDynamicsStream, trialTimeSeries, calculateIndex, meanDistanceMatrix, nearestNeighbors, useLocalDimensions, minReturnTime);

        waitHandle.iterate(1);
    end
    close(waitHandle);

    plotIndices = calcIndicies;
    bestLocalDistances = localDistances';

    smoothedNeighbors = filterData(bestNeighborCounts, 3);
        
    %% Use in degree vs out degree to find bad points

    figure(FIG_BASE + 2);
    clf;
    h(1) = subplot(2,1,1);
    imagesc(bestLocalProjections);
    title('Local projections (findBestAxes)');
    ylabel('Projection dim');
    h(2) = subplot(2,1,2);
    imagesc(thisDynamicsStream');
    title('Dynamics stream (neurons x time)');
    xlabel('Time index');
    ylabel('Neuron index');
    linkaxes(h);
    colormap(parula(256));

    test = bestLocalDistances ./ sum(bestLocalDistances,2);
    % clf
    % plot(max(test,[],2)' ./ sum(test));

    inCounts = [];
    outCounts = [];
    for i = 1:size(bestLocalDistances,1)
        thisProbabilities = bestLocalDistances(i,:) > exp(-1);
        outCounts(i) = max(bwlabel(thisProbabilities));

        thisProbabilities = bestLocalDistances(:,i) > exp(-1);
        inCounts(i) = max(bwlabel(thisProbabilities));
    end

    allStateValidities = inCounts ./ outCounts;

    figure(FIG_BASE + 3);
    clf;
    plot(allStateValidities);
    title('In/Out degree ratio (state validities)');
    xlabel('State index');
    ylabel('In/Out ratio');

%%

    IN_OUT_RATIO_CUTOFF = 0;

    diffusionMapIndicies = 1:size(bestLocalDistances,1);

    MIN_TRAJECTORY_SIZE = 5;
    TRAJECTORY_SIZE = 2;
    x = -3*TRAJECTORY_SIZE:3*TRAJECTORY_SIZE;
    gaussian = exp(-x.^2/(2*TRAJECTORY_SIZE^2));
    gaussian = gaussian / sum(gaussian);
    kernel = diag(gaussian);

    smoothedLocalDistances = bestLocalDistances;
    smoothedLocalDistances = min(smoothedLocalDistances, smoothedLocalDistances');

    trimmedMatrix = smoothedLocalDistances;

    terminalIDs = [];
    currentID = 0;
    for i = 1:length(trialTimeSeries)
        currentID = currentID + size(trialTimeSeries{i},2);
        terminalIDs = [terminalIDs, currentID];
    end

    allTerminalIDs = terminalIDs;

    trimmedMatrix(terminalIDs,:) = 0;
    trimmedMatrix(:,terminalIDs) = 0;   

    trimmedMatrix(trimmedMatrix < exp(-2)) = 0;
    goodIndices = find(sum(trimmedMatrix,2) ~= 0);
    badIndices = find(sum(trimmedMatrix,2) == 0);
    for i = 1:length(badIndices)
        trimmedMatrix(badIndices(i),badIndices(i)) = 1;
    end
    trimmedMatrix(goodIndices,:) = trimmedMatrix(goodIndices,:) ./ nansum(trimmedMatrix(goodIndices,:), 2);
    exponentiatedMatrix = trimmedMatrix;
    exponentCount = 1;
    lastPopulation = -1;
    while max(sum(exponentiatedMatrix(goodIndices,:) > 0, 2)) < length(goodIndices) * repopulateDensity && max(sum(exponentiatedMatrix(goodIndices,:) > 0, 2)) ~= lastPopulation
        exponentCount = exponentCount + 1;

        lastPopulation = max(sum(exponentiatedMatrix(goodIndices,:) > 0, 2));

        exponentiatedMatrix = exponentiatedMatrix * trimmedMatrix;
    end

    trimmedMatrix = exponentiatedMatrix;

    allStarts = [1 terminalIDs+1];
    allStarts(end) = [];

    thisStart = 1;
    for i = 1:length(terminalIDs)
        if i > 1
            thisStart = terminalIDs(i-1) + 1;
        end
        thisEnd = terminalIDs(i);

        trimmedMatrix(thisStart:thisEnd-1, thisStart:thisEnd-1) = trimmedMatrix(thisStart+1:thisEnd, thisStart:thisEnd-1);
        trimmedMatrix(thisEnd-1,:) = 0;
        trimmedMatrix(thisEnd-1,allStarts) = 1;
    end

    trimmedMatrix(terminalIDs,:) = [];
    trimmedMatrix(:,terminalIDs) = [];   

    exponentiatedMatrix = trimmedMatrix;

    terminalIDs = allTerminalIDs - (1:length(allTerminalIDs));
    allTerminalIDs = [allTerminalIDs allTerminalIDs-1];

    asymmetricProbabilities = exponentiatedMatrix;%(2:end,1:end-1);
    asymmetricProbabilities(terminalIDs,:) = [];
    asymmetricProbabilities(:,terminalIDs) = [];
    asymmetricProbabilities = asymmetricProbabilities ./ sum(asymmetricProbabilities, 2);

    [steadyState,~] = eigs(asymmetricProbabilities',1,'largestabs','Tolerance',1e-16,'MaxIterations',1000);
    steadyState = steadyState / sum(steadyState);
    
    if sum(steadyState(:) == 0) > 0
        steadyState = sum(asymmetricProbabilities,1); %If eigenvalue estimation is unstable fall back to incoming degree
    end

    stateDensities = steadyState * steadyState';
    stateDensities(stateDensities < 0) = 0; % Rounding errors can give <0 and cause complex numbers
    

    
    stateDensities = sqrt(stateDensities);

    asymmetricProbabilities = asymmetricProbabilities ./ stateDensities;

    asymmetricProbabilities = asymmetricProbabilities ./ sum(asymmetricProbabilities,2);

    allBadIndices = [];

    badIndicies =  find(inCounts ./ outCounts < IN_OUT_RATIO_CUTOFF);        
    badIndicies(badIndicies > size(asymmetricProbabilities,1)) = [];

    asymmetricProbabilities(badIndicies,:) = 0;
    asymmetricProbabilities(:,badIndicies) = 0;

    allBadIndices = badIndicies;

    reducedProbabilities = sum(asymmetricProbabilities,2);        
    badIndicies = find(reducedProbabilities == 0);

    while length(badIndicies) > length(allBadIndices)
        asymmetricProbabilities(badIndicies,:) = 0;
        asymmetricProbabilities(:,badIndicies) = 0;

        allBadIndices = unique([allBadIndices, badIndicies']);

        reducedProbabilities = sum(asymmetricProbabilities,2);        
        badIndicies = find(reducedProbabilities == 0);
    end

    asymmetricProbabilities(allBadIndices,:) = [];
    asymmetricProbabilities(:,allBadIndices) = [];

    finalIndicies = 1:size(thisDynamicsStream,1);
    finalIndicies(allBadIndices) = [];

    finalDynamicsStream = thisDynamicsStream(1:end,:);
    finalDynamicsStream(allTerminalIDs,:) = [];    

    stateValidities = allStateValidities;
    stateValidities(allBadIndices) = [];    

%         diffusionMapIndicies(end) = [];
    diffusionMapIndicies(allBadIndices) = [];

    stateHasNext = diff(diffusionMapIndicies) == 1;

    jumpIndicies = find(stateHasNext == 0);

    figure(FIG_BASE + 4);
    clf;
    hold on;
    if size(pcaBasis,2) > 2
        h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
    else
        h = plot(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), 'LineWidth', 0.5);
    end
    h.Color(4) = 0.2;
    title('Final dynamics stream (PCA)');
    if size(pcaBasis,2) > 2
        scatter3(finalDynamicsStream(jumpIndicies,:)*pcaBasis(:,1), finalDynamicsStream(jumpIndicies,:)*pcaBasis(:,2), finalDynamicsStream(jumpIndicies,:)*pcaBasis(:,3), 32, 'kx');
        xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
    else
        scatter(finalDynamicsStream(jumpIndicies,:)*pcaBasis(:,1), finalDynamicsStream(jumpIndicies,:)*pcaBasis(:,2), 32, 'kx');
        xlabel('PC1'); ylabel('PC2');
    end
