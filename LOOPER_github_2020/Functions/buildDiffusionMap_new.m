
%% get variables

% usePoints = [1000:6000];%; 10000:15000; 20000:25000];%; 30000:35000; 40000:45000];
% useTrialIDs = 1;%[1,4; 2,5; 3,6];%; 4,7; 5,8];
% savePrefixs = {};  
% for dataSetIndex = 1:size(usePoints,1)
%     savePrefixs{dataSetIndex} = ['systemCTRNN_NewSigma' num2str(dataSetIndex)];
% end
% % usePoints = [1000:6000];
% timeSkip = 1;
% 
% for dataSetIndex = 1%:size(usePoints,1)

%     inputStream = inputs(usePoints(dataSetIndex,:));
%     targetStream = targets(usePoints(dataSetIndex,:),:);
%     % outputStream = outputs(usePoints);
%     % networkOutputStream = testY(usePoints);
%     rawDynamicsStream = dynamics(usePoints(dataSetIndex,:),:,:);
%     rawDynamicsStream = permute(rawDynamicsStream, [1, 3, 2]);
%     % allDynamicsStream = [allDynamicsStream, normrnd(0,10,[size(allDynamicsStream,1), 100, 10])];
%     totalTrials = 1;%size(rawDynamicsStream,3);
%     
%     noiseData = [];
%     noiseData(1) = 0;
%     for j = 2:100000
%         noiseData(j) = noiseData(j-1) * 0.95 + normrnd(0, 10, 1, 1);
%     end
%     
%     augmentedDynamics = [];
%     for i = 1:totalTrials
% %         noiseLocations = randsample(size(rawDynamicsStream, 1) + 1:length(noiseData), 2000);
% %         noiseNeurons = [];
% %         for j = 1:length(noiseLocations)
% %             noiseNeurons(:,j) = noiseData(noiseLocations(j)+(1-size(rawDynamicsStream, 1):0));
% %         end
% %         noiseNeurons = repmat(noiseNeurons, [1 100]);
% %         gaussianNeurons = normrnd(0, 100, [size(rawDynamicsStream,1), 1000]);
% %         gaussianNeurons = repmat(gaussianNeurons, [1 30]);
%         
%         noiseNeurons = [];
%         gaussianNeurons = [];
%         
%         augmentedDynamics(:,:,i) = zscore([rawDynamicsStream(:,:,i), gaussianNeurons, noiseNeurons]);
%     end
%     
%     allDynamics = [];
%     for i = 1:totalTrials
%         allDynamics = [allDynamics; augmentedDynamics(:,:,i)];
%     end
%     [fullBasis, ~] = pca(allDynamics);
% 
%     embeddedTime = 0;
%     allDynamicsStream = [];
%     for i = 1:totalTrials
%         tempDynamics = squeeze(augmentedDynamics(:,:,i));% * fullBasis;
%         
% %         tempDynamics = filterData(tempDynamics, 1);
% 
% %         [tempDynamics, embeddedTime] = delayEmbed(tempDynamics, 2, 1, 0);
% 
%         thisIndices = 1:size(tempDynamics,1);
%     %     tempDynamics = [tempDynamics, targetStream(thisIndices+embeddedTime,2)];
%         tempDynamics = [tempDynamics, inputStream(thisIndices+embeddedTime,:)];
% 
%         allDynamicsStream(:,:,i) = tempDynamics;
%     end
% 
% 
% 
% 
    
% 
% 
%     thisDynamicsStream = allDynamicsStream(:,:,1);
%     
%     goodDimensions = 1:size(thisDynamicsStream,2);
% %     goodDimensions = 1:100;
% 
%     [pcaBasis, pcaOutputs] = pca(thisDynamicsStream(:,goodDimensions), 'NumComponents', 10);
% 
%     plotIndices = 1:size(thisDynamicsStream,1);
%     % plotIndices = size(allDynamicsStream,1)-1000:size(allDynamicsStream,1);
% 
%     % stateIDs = targetStream(plotIndices,2) + (targetStream(plotIndices,1)>exp(-1))*4;
%     % stateIDs = stateIDs + 1;
% 
%     targetStream = targetStream(plotIndices+embeddedTime,:);
%     inputStream = inputStream(plotIndices+embeddedTime,:);
%     times = 1:length(inputStream);
% 
%     stateIDs = [];
%     currentState = 1;
%     for i = 1:size(targetStream,1)
%         if currentState == 1
%             if inputStream(i,1) == 2
%                 currentState = 2;
%             end
%         else
%             if inputStream(i,1) == 1
%                 currentState = 1;
%             end
%         end
% 
%         stateIDs(i) = currentState;
%     end
%     % stateIDs = stateIDs' + (targetStream(plotIndices,1)>exp(-1))*2;
%     stateIDs = stateIDs' + (targetStream(plotIndices,1))*2;
%     % stateIDs = (targetStream(plotIndices,1)>exp(-1))+1;
% 
%     PCMode = 2;
% 
%     figure(1);
%     clf;
%     colors = lines(8);
%     plotColor3(thisDynamicsStream(plotIndices,goodDimensions)*pcaBasis(:,1),thisDynamicsStream(plotIndices,goodDimensions)*pcaBasis(:,2),thisDynamicsStream(plotIndices,goodDimensions)*pcaBasis(:,3), stateIDs, colors)
%     xlabel('DPC1');
%     ylabel('DPC2');
%     zlabel('DPC3');
% 
% %     figure(2);
% %     clf;
% %     colors = lines(8);
% %     plotColor(plotIndices, thisDynamicsStream(plotIndices,goodDimensions)*pcaBasis(:,PCMode), stateIDs, colors)
% %     hold on;
% %     plot(inputStream(plotIndices), 'k')
% %     ylabel('PC / Input');
% %     xlabel('Timesteps');
%     
%     figure(3);
%     clf;
%     imagesc(thisDynamicsStream);
    
    
    %%
%     for trialIDIndex = 1%:length(useTrialIDs(dataSetIndex,:))
        %%
%         trialID = useTrialIDs(dataSetIndex, trialIDIndex);        
        thisDynamicsStream = timeSeriesData';
        
        [pcaBasis, pcaOutputs] = pca(thisDynamicsStream, 'NumComponents', 3);
        
%         derivatives = [diff(thisDynamicsStream); zeros(size(thisDynamicsStream(1,:)))];
        
%         thisDynamicsStream = [thisDynamicsStream derivatives];
%         meanDistanceMatrix = distanceMatrix(:,:,trialID);
        meanDistanceMatrix = pdist2(thisDynamicsStream, thisDynamicsStream);
%         meanDistanceMatrix = pdist2(thisDynamicsStream, thisDynamicsStream);

        bestSigmaValues = [];
        bestNeighborCounts = [];
        bestLocalProjections = [];
%         bestLocalDistances = [];
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
        
%         smoothSigma = 1;
%         x = -ceil(3*smoothSigma):ceil(3*smoothSigma);
%         kernel = -x.*exp(-x.^2/(2*smoothSigma^2))/(smoothSigma^3*sqrt(2*pi));
% 
%         dynamicsDiff = [];
%         for j = 1:length(trialTimeSeries)
%             thisTrial = trialTimeSeries{j}';
%             dynamicsDiff = [dynamicsDiff; filterData(thisTrial, 0, kernel, 1, 0)];
%             dynamicsDiff(end,:) = nan(size(dynamicsDiff(1,:)));
%         end

        waitHandle = parfor_progressbar(size(thisDynamicsStream,1), 'Calculating nearest neighbors');        
%         while min(bestDistances) < 5
        parfor calculateIndex = calcIndicies
%             calculateIndex = randsample(calcIndicies, 1, true, 1./bestDistances);
            
%             [bestLocalProjections(:,calculateIndex), bestSigmaValues(calculateIndex), bestDiffSigmaValues(calculateIndex), localDistances(:,calculateIndex), bestNeighborCounts(calculateIndex), distances] = findBestAxes(thisDynamicsStream, dynamicsDiff, trialTimeSeries, calculateIndex, meanDistanceMatrix, nearestNeighbors, useLocalDimensions);
            [bestLocalProjections(:,calculateIndex), bestSigmaValues(calculateIndex), bestDiffSigmaValues(calculateIndex), localDistances(:,calculateIndex), bestNeighborCounts(calculateIndex), distances] = findBestAxes(thisDynamicsStream, trialTimeSeries, calculateIndex, meanDistanceMatrix, nearestNeighbors, useLocalDimensions);
        
%             replacedIndices = (distances > bestDistances);
%             
%             bestDistances(replacedIndices) = distances(replacedIndices);
%             bestNeighborhood(replacedIndices) = calculateIndex;
            
            waitHandle.iterate(1);
        end
        close(waitHandle);
        
        plotIndices = calcIndicies;
        bestLocalDistances = localDistances';
        
        smoothedNeighbors = filterData(bestNeighborCounts, 3);
        
%         %% Rebuild distances
%         
%         FILTER_AMOUNT = 1;
%         
%         smoothedLocalProjections = [];
%         for i = 1:size(bestLocalProjections,1)
%             smoothedLocalProjections(i,:) = filterData(bestLocalProjections(i,:), FILTER_AMOUNT, [], false, 0);
%         end
%         smoothedLocalProjections = smoothedLocalProjections ./ sqrt(sum(smoothedLocalProjections.^2));
%         
%         smoothedSigmas = filterData(bestSigmaValues, FILTER_AMOUNT);
%         smoothedSigmas = smoothedSigmas;
%         
%         
%         figure(4);
%         clf;
%         h(1) = subplot(2,1,1);
%         imagesc(smoothedLocalProjections);
%         h(2) = subplot(2,1,2);
%         imagesc(thisDynamicsStream');
%         linkaxes(h);
%         colormap(parula(256));
%         
%         %%
%         
%         smoothedNeighbors = [];
%         localDistances = [];
%         waitHandle = parfor_progressbar(length(smoothedSigmas), 'Building diffusion map');
%         parfor i = 1:length(smoothedSigmas)
% %         for i = 25
%             localStream = thisDynamicsStream ./ bestLocalProjections(:,i)';
% %             localStream = thisDynamicsStream ./ smoothedLocalProjections(:,i)';
%             thisValues = pdist2(localStream(i,:), localStream);
%                 
%             smoothSigma = 1;
%             x = -ceil(3*smoothSigma):ceil(3*smoothSigma);
%             kernel = -x.*exp(-x.^2/(2*smoothSigma^2))/(smoothSigma^3*sqrt(2*pi));
% 
%             currentDiff = filterData(localStream, 0, kernel, 1, 0);
%             thisDiffDistances = pdist2(currentDiff(showIndex,:), currentDiff, 'cosine');
%             sortedDiffDistances = sort(thisDiffDistances);
%             diffCurvature = filterData(filterData(sortedDiffDistances, 0, kernel, 1, 0), 50);
%             [~, bestDiffSigmaIndex] = min(diffCurvature);
%             bestDiffSigma = sortedDiffDistances(ceil(bestDiffSigmaIndex/2));
%             guassianDiff = exp(-thisDiffDistances.^2/(2*bestDiffSigma^2));
%             guassianDiff(guassianDiff < 0.1) = 0.1;
% 
%             thisValues = thisValues ./ guassianDiff;
%             
%             [peaks, peakIDs] = findpeaks(-thisValues);
%             peaks = -peaks;
% 
%             [sortedPeaks, sortIDs] = sort(peaks);
% 
%             sortedPeakTimes = peakIDs(sortIDs);
% 
%             j = 1;
%             while j < length(sortedPeakTimes)
%                 tempPeakTimes = sortedPeakTimes(j+1:end);
%                 repeatIndices = find(tempPeakTimes > sortedPeakTimes(j) - MIN_RETURN_TIME & tempPeakTimes < sortedPeakTimes(j) + MIN_RETURN_TIME);
% 
%                 sortedPeakTimes(repeatIndices + j) = [];
% 
% 
%                 j = j + 1;
%             end
%             
%             allValues = exp(-thisValues.^2/(2*(bestSigmaValues(i))^2));
%             localDistances(i,:) = allValues;
%             smoothedNeighbors(i) = sum(allValues(sortedPeakTimes) > exp(-1));
%             
%             waitHandle.iterate(1);
%         end
%         close(waitHandle);
%         
        
        
        
        %%
        
%         figure(4);
%         clf;
%         hold on;
%         colors = lines(8);
%         h = plot3(plotStream*pcaBasis(:,1),plotStream*pcaBasis(:,2),plotStream*pcaBasis(:,3));
%         h.Color(4) = 0.2;
%         scatter3(plotStream(plotIndices,goodDimensions)*pcaBasis(:,1),plotStream(plotIndices,goodDimensions)*pcaBasis(:,2),plotStream(plotIndices,goodDimensions)*pcaBasis(:,3), 32, bestSigmaValues(plotIndices))
%         xlabel('DPC1');
%         ylabel('DPC2');
%         zlabel('DPC3');
%         colormap(jet(256));

        %% Get detailed balance decomp
        
%         totalValues = sum(localDistances,2);
%         tempProbabilities = localDistances ./ totalValues;
%         tempProbabilities = tempProbabilities^4;
%         
%         [steadyState, ~] = eigs(tempProbabilities', 1);
%         fluxMatrix = diag(steadyState) * tempProbabilities;
%         symmetricMatrix = (fluxMatrix+fluxMatrix.')/2;
%         antisymmetricMatrix = (fluxMatrix-fluxMatrix.')/2;
%         symmetricMatrix = symmetricMatrix ./ steadyState;
%         antisymmetricMatrix = antisymmetricMatrix ./ steadyState;
% 
%         diffusedProbabilities = max(0, symmetricMatrix^4 + antisymmetricMatrix);
%         diffusedProbabilities = max(0, symmetricMatrix + antisymmetricMatrix);
%         diffusedProbabilities = diffusedProbabilities ./ sum(diffusedProbabilities,2);
%         
%         diffusedDistances = diffusedProbabilities .* totalValues;
        
        %% Use in degree vs out degree to find bad points
        
        figure(2);
        clf;
        h(1) = subplot(2,1,1);
        imagesc(bestLocalProjections);
        h(2) = subplot(2,1,2);
        imagesc(thisDynamicsStream');
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

        clf
        plot(allStateValidities);
        
%%

%         IN_OUT_RATIO_CUTOFF = 0.6;
        IN_OUT_RATIO_CUTOFF = 0;
        
        diffusionMapIndicies = 1:size(bestLocalDistances,1);

        MIN_TRAJECTORY_SIZE = 5;
        TRAJECTORY_SIZE = 2;
        x = -3*TRAJECTORY_SIZE:3*TRAJECTORY_SIZE;
        gaussian = exp(-x.^2/(2*TRAJECTORY_SIZE^2));
        gaussian = gaussian / sum(gaussian);
        kernel = diag(gaussian);
        
%         smoothedLocalDistances = convn(bestLocalDistances, kernel, 'same');
        smoothedLocalDistances = bestLocalDistances;
%         smoothedLocalDistances = (smoothedLocalDistances - exp(-1)) / (1 - exp(-1));
        smoothedLocalDistances = min(smoothedLocalDistances, smoothedLocalDistances');
%         smoothedLocalDistances = (smoothedLocalDistances + smoothedLocalDistances')/2;
%         smoothedLocalDistances = sqrt(smoothedLocalDistances .* smoothedLocalDistances');
%         smoothedLocalDistances = (smoothedLocalDistances - exp(-1)) / (1 - exp(-1));
%         smoothedLocalDistances(smoothedLocalDistances < 0) = 0;
%         smoothedLocalDistances(smoothedLocalDistances < exp(-1)) = 0;
%         smoothedLocalDistances = (smoothedLocalDistances).^2;
%         smoothedLocalDistances = min(smoothedLocalDistances, smoothedLocalDistances');
%         smoothedLocalDistances = sqrt(smoothedLocalDistances * smoothedLocalDistances');
%         smoothedLocalDistances = smoothedLocalDistances ./ sum(smoothedLocalDistances,2);

%         smoothedLocalDistances = (smoothedLocalDistances - exp(-1)) / (1 - exp(-1));

%         finalNeighbors = [];
%         for i = 1:size(smoothedProbabilities)
%             thisProbabilities = smoothedLocalDistances(i,:) > exp(-1);
%             finalNeighbors(i) = max(bwlabel(thisProbabilities));
%         end
%         finalNeighbors = filterData(finalNeighbors, 3);

%         stateDensities = sum(smoothedLocalDistances);
%         
%         stateDensities = (stateDensities' * stateDensities);
%         smoothedLocalDistances = smoothedLocalDistances ./ stateDensities;
% 
%         totalValue = sum(smoothedLocalDistances,2);
%         smoothedProbabilities = smoothedLocalDistances ./ totalValue;
%         maxProbabilities = max(smoothedProbabilities,[],2);
%         badIndicies = maxProbabilities > 0.25;

%         mixedProbabilities = smoothedProbabilities^2;
%         smoothedProbabilities(mixedProbabilities < exp(-1))

%         MAX_PROBABILITY = 0.2;
        
        trimmedMatrix = smoothedLocalDistances;
        terminalIDs = find(isnan(trimmedMatrix(1,:)));
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
%         exponentiatedMatrix = trimmedMatrix;
%         exponentCount = 0;
% %         while min(sum(exponentiatedMatrix(goodIndices,:) > 0, 2)) < size(exponentiatedMatrix,1) / 4 && exponentCount < 10
% %             exponentCount = exponentCount + 1;
% %             exponentiatedMatrix = exponentiatedMatrix * trimmedMatrix;
% %         end
%         
%         trimmedMatrix = exponentiatedMatrix;
        
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
        
        
        
% %         trimmedMatrix(end,:) = 0;
% %         trimmedMatrix(:,end) = 0;
%         
% %         trimmedMatrix(trimmedMatrix < exp(-2)) = 0;
%         goodIndices = find(sum(trimmedMatrix,2) ~= 0);
%         badIndices = find(sum(trimmedMatrix,2) == 0);
%         for i = 1:length(badIndices)
%             trimmedMatrix(badIndices(i),badIndices(i)) = 1;
%         end
%         trimmedMatrix(goodIndices,:) = trimmedMatrix(goodIndices,:) ./ nansum(trimmedMatrix(goodIndices,:), 2);
%         exponentiatedMatrix = trimmedMatrix;
%         exponentCount = 1;
%         while min(sum(exponentiatedMatrix(goodIndices,:) > 0, 2)) < size(exponentiatedMatrix,1) / 4 && exponentCount < 10
%             exponentCount = exponentCount + 1;
%             exponentiatedMatrix = exponentiatedMatrix * trimmedMatrix;
%         end
%         
% %         for i = 1:length(badIndices)
% %             exponentiatedMatrix(badIndices(i),badIndices(i)) = 0;
% %         end

        terminalIDs = allTerminalIDs - (1:length(allTerminalIDs));
        allTerminalIDs = [allTerminalIDs allTerminalIDs-1];
        
        asymmetricProbabilities = trimmedMatrix;%(2:end,1:end-1);
        asymmetricProbabilities(terminalIDs,:) = [];
        asymmetricProbabilities(:,terminalIDs) = [];
%         badIndices = find(sum(asymmetricProbabilities,2) == 0);
%         for i = 1:length(asymmetricProbabilities)
%             asymmetricProbabilities(badIndices,badIndices) = 1;
%         end

        asymmetricMarkov = asymmetricProbabilities ./ sum(asymmetricProbabilities, 2);
        
        [steadyState,~] = eigs(asymmetricMarkov',1,'largestabs','Tolerance',1e-16,'MaxIterations',1000);
        steadyState = steadyState / sum(steadyState);

        if var(steadyState) == 0
            steadyState = sum(asymmetricProbabilities^100,1) / size(asymmetricProbabilities,1);
        end
        
        stateDensities = sqrt(steadyState * steadyState');
        asymmetricProbabilities = asymmetricProbabilities ./ stateDensities;

        asymmetricProbabilities = asymmetricProbabilities ./ sum(asymmetricProbabilities,2);
        
%         [steadyState,~] = eigs(asymmetricProbabilities',1);
%         steadyState = steadyState / sum(steadyState);         
%         figure(3);
%         clf;
%         title('Laplace Beltrami');
%         hold on;
%         h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
%         h.Color(4) = 0.2;
%         scatter3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 32, steadyState);
%         colormap(jet(256));
        
        allBadIndices = [];
        
%         maxProbabilities = max(smoothedProbabilities,[],2);
%         badIndiciesRow = maxProbabilities > MAX_PROBABILITY;
%         maxProbabilities = max(smoothedProbabilities,[],1);
%         badIndiciesColumn = maxProbabilities > MAX_PROBABILITY;
%         badIndicies = unique([find(badIndiciesRow); find(badIndiciesColumn)']);
%         badIndicies = finalNeighbors <= 5;
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
        
%         badIndiciesLogical = 0 * (1:size(asymmetricProbabilities,1));
%         badIndiciesLogical(allBadIndices) = 1;
%         
% %         badIndiciesLogical = convn(badIndiciesLogical, ones(1,MIN_TRAJECTORY_SIZE), 'same');
% %         badIndiciesLogical = badIndiciesLogical > 0;
%         
%         badIndiciesLogical = imclose(badIndiciesLogical, ones(1,MIN_TRAJECTORY_SIZE));
% %         badIndiciesLogical = imopen(badIndiciesLogical, ones(1,MIN_TRAJECTORY_SIZE));
%         
%         badIndicies = find(badIndiciesLogical);
%         
%         asymmetricProbabilities(badIndicies,:) = 0;
%         asymmetricProbabilities(:,badIndicies) = 0;
%         
%         reducedProbabilities = sum(asymmetricProbabilities,2);        
%         badIndicies = find(reducedProbabilities == 0);
%         
%         while length(badIndicies) > length(allBadIndices)
%             asymmetricProbabilities(badIndicies,:) = 0;
%             asymmetricProbabilities(:,badIndicies) = 0;
%             
%             allBadIndices = unique([allBadIndices, badIndicies']);
%             
%             reducedProbabilities = sum(asymmetricProbabilities,2);        
%             badIndicies = find(reducedProbabilities == 0);
%         end
        
        asymmetricProbabilities(allBadIndices,:) = [];
        asymmetricProbabilities(:,allBadIndices) = [];
        
%         localWeights = sum(asymmetricProbabilities);
%         asymmetricProbabilities = asymmetricProbabilities ./ sqrt(localWeights'*localWeights);
%         asymmetricProbabilities = asymmetricProbabilities ./ sum(asymmetricProbabilities, 2);


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
        figure(3);
        clf;
        hold on;
        h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
        h.Color(4) = 0.2;
        h = plot3(thisDynamicsStream*pcaBasis(:,1), thisDynamicsStream*pcaBasis(:,2), thisDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
        h.Color(4) = 0.2;
        scatter3(finalDynamicsStream(jumpIndicies,:)*pcaBasis(:,1), finalDynamicsStream(jumpIndicies,:)*pcaBasis(:,2), finalDynamicsStream(jumpIndicies,:)*pcaBasis(:,3), 32, 'kx');

        
%         asymmetricProbabilities = asymmetricProbabilities ./ sum(asymmetricProbabilities,2);
%         maxProbabilities = max(asymmetricProbabilities,[],2);
%         
%         while max(maxProbabilities) > 0.60
%             badIndiciesRow = find(maxProbabilities > 0.60);
%             asymmetricProbabilities(badIndiciesRow,:) = [];
%             asymmetricProbabilities(:,badIndiciesRow) = [];
% 
%             asymmetricProbabilities = asymmetricProbabilities ./ sum(asymmetricProbabilities,2);
%             maxProbabilities = max(asymmetricProbabilities,[],2);
% 
%             finalDynamicsStream(badIndiciesRow,:) = [];    
%         end
        
        
%         nearTrajectoryCounts = [];
%         for i = 1:length(smoothedSigmas)
%             values = smoothedLocalDistances(i,:);
%             showIndicies = find(values > exp(-1));
%             
%             j = 1;
%             while j < length(showIndicies)
%                 tempPeakTimes = showIndicies(j+1:end);
%                 repeatIndices = find(tempPeakTimes > showIndicies(j) - MIN_RETURN_TIME & tempPeakTimes < showIndicies(j) + MIN_RETURN_TIME);
% 
%                 showIndicies(repeatIndices + j) = [];
% 
%                 j = j + 1;
%             end
%         
%             nearTrajectoryCounts(i) = length(showIndicies);
%         end
        
%         badIndicies = [];
%         badIndicies = nearTrajectoryCounts < 5;
        
%         smoothedLocalDistances(badIndicies,badIndicies) = nan;

%         validConnections = smoothedLocalDistances;% * smoothedProbabilities';
%         for i = 1:size(validConnections,1)
%             validConnections(i,:) = validConnections(i,:) .* validConnections(:,i)';
%             validConnections(i,validConnections(i,:) < exp(-2)) = 0;
%             minValue = min(validConnections(i,:));
%             maxValue = max(validConnections(i,:));
%             validConnections(i,:) = (validConnections(i,:) - minValue) / (maxValue - minValue);
%             
%             validConnections(i,:) = validConnections(i,:) > 0;
%         end        
%         validConnections = min(validConnections, validConnections');
%         
%         smoothedProbabilities(~validConnections) = 0;
%         
%         figure(1);
%         clf;
%         imagesc(smoothedLocalDistances(calcIndicies,:))
        
%         values = maxProbabilities;
%         values(badIndicies) = nan;
%         
%         plotIndices = 1:size(plotStream,1);
        
        
%         figure(4);
%         clf;
%         hold on;
%         colors = lines(8);
%         h = plot3(plotStream*pcaBasis(:,1),plotStream*pcaBasis(:,2),plotStream*pcaBasis(:,3));
%         h.Color(4) = 0.2;
%         scatter3(plotStream(plotIndices,goodDimensions)*pcaBasis(:,1),plotStream(plotIndices,goodDimensions)*pcaBasis(:,2),plotStream(plotIndices,goodDimensions)*pcaBasis(:,3), 32, values(plotIndices))
%         xlabel('DPC1');
%         ylabel('DPC2');
%         zlabel('DPC3');
%         colormap(jet(256));
%         caxis([0 30]);
        
        %%
        
%         [eigenvector, eigvalue] = eigs(smoothedProbabilities', 1);
%         %%
%         %[caz,cel] = view;
%         
%         displayIndicies = 1:size(thisDynamicsStream,1);
% %         displayIndicies = 350:360;
%         
% %         for i = 1:5001
%         for i = 3461
%             displayIndex = i;
% 
%     %         tests = smoothedLocalDistances;
% %             test = smoothedProbabilities * smoothedProbabilities';
% 
% %             values = validConnections(displayIndex,:);
%             showIndicies = find(smoothedProbabilities(displayIndex,:) > 0.01);
%             values = smoothedProbabilities(displayIndex,:);
%             
% %             showIndicies = 1:size(thisDynamicsStream,1);
% %             values = max(test,[],2)' ./ sum(test) .* finalNeighbors.^2;
% 
%             figure(5);
%             clf;
%             hold on;
%             colors = lines(8);
%             h = plot3(thisDynamicsStream(displayIndicies,:)*pcaBasis(:,1),thisDynamicsStream(displayIndicies,:)*pcaBasis(:,2),thisDynamicsStream(displayIndicies,:)*pcaBasis(:,3));
%             h.Color(4) = 0.2;
%             scatter3(thisDynamicsStream(showIndicies,:)*pcaBasis(:,1),thisDynamicsStream(showIndicies,:)*pcaBasis(:,2),thisDynamicsStream(showIndicies,:)*pcaBasis(:,3), 32, values(showIndicies))
%             xlabel('DPC1');
%             ylabel('DPC2');
%             zlabel('DPC3');
%             colormap(jet(256));
%             title(num2str(i));
%             if exist('caz')
%                 view(caz,cel)
%             end
%             
%             drawnow();
%             
% %             pause(0.2 / length(showIndicies));
%         end
% %     end
%     
% % end
