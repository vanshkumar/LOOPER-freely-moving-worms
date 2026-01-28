function [bestLocalProjections, bestSigmaValues, bestDiffSigmaValues, localDistances, bestNeighborCounts, distances] = findBestAxes(dynamicsStream, trialStream, calcIndex, meanDistanceMatrix, nearestNeighbors, useLocalDimensions, minReturnTime, verbose)


    bestSigmaValues = [];
    bestNeighborCounts = [];
    bestLocalProjections = [];
    distances = [];
    bestDiffSigmaValues = [];
    localDistances = [];

    %% Near trajectories
    
    if ~exist('meanDistanceMatrix', 'var') || isempty(meanDistanceMatrix)
        meanDistanceMatrix = pdist2(dynamicsStream, dynamicsStream);
    end
    if ~exist('verbose', 'var') || isempty(verbose)
        verbose = 0;
    end
    if ~exist('nearestNeighbors', 'var') || isempty(nearestNeighbors)
        nearestNeighbors = 12;
    end
    if ~exist('useLocalDimensions', 'var') || isempty(useLocalDimensions)
        useLocalDimensions = 1;
    end
    if ~exist('minReturnTime', 'var') || isempty(minReturnTime)
        minReturnTime = 10;
    end

tic
    
    test = meanDistanceMatrix(1:end,1:end);
    totalValues = test(calcIndex,:);

    [peaks, peakIDs] = findpeaks(-totalValues);
    peaks = -peaks;

    [sortedPeaks, sortIDs] = sort(peaks);

    sortedPeakTimes = peakIDs(sortIDs);
    
    trialEnds = [];
    lastEnd = 0;
    for i = 1:length(trialStream)
        lastEnd = lastEnd + size(trialStream{i},2);
        trialEnds(i) = lastEnd;
    end    
    trialStarts = [1 trialEnds+1];

    i = 1;
    while i < length(sortedPeakTimes)
        tempPeakTimes = sortedPeakTimes(i+1:end);
        repeatIndices = find(tempPeakTimes > sortedPeakTimes(i) - minReturnTime & tempPeakTimes < sortedPeakTimes(i) + minReturnTime);

        sortedPeakTimes(repeatIndices + i) = [];


        i = i + 1;
    end

    %%

    MIN_NEIGHBOR_COUNT = 0;
    BEST_TRAJECTORY_COUNT = 1;

    currentPeakCount = nearestNeighbors;

    maxDistances = [];
    bestSigmas = [];
    bestProjections = [];
    bestDistances = [];
    peakIndices = [];
    peakValues = [];
    bestValues = [];
    worstDistances = [];
    minDistances = [];
    neighborhoodCounts = [];
    diffSigmas = [];
    gaussianDistances = [];
    
    for i = 1
            currentPeaks = peakIDs(sortIDs(1:currentPeakCount));
        stdPeaks = [currentPeaks-1 currentPeaks+1];
        stdPeaks(stdPeaks < 1) = [];
        stdPeaks(stdPeaks > size(dynamicsStream,1)) = [];
        stdPeaks(ismember(stdPeaks, trialStarts)) = [];
        stdPeaks(ismember(stdPeaks, trialEnds)) = [];
        stdPeaks = [currentPeaks, stdPeaks];
        
        if useLocalDimensions
            currentSTD = std(dynamicsStream(stdPeaks,:));
        else
            currentSTD = ones([1, size(dynamicsStream,2)]);
        end
        currentSTD = currentSTD / sqrt(sum(currentSTD.^2));
        currentSTD(currentSTD < 1e-4) = 1e-4;
        currentStream = dynamicsStream ./ currentSTD;
        
        thisValues = pdist2(currentStream(currentPeaks,:), currentStream);

        smoothSigma = 1;
        x = -ceil(3*smoothSigma):ceil(3*smoothSigma);
        kernel = -x.*exp(-x.^2/(2*smoothSigma^2))/(smoothSigma^3*sqrt(2*pi));

        currentDiff = [];
        for j = 1:length(trialStream)
            thisTrial = trialStream{j}' ./ currentSTD;
            currentDiff = [currentDiff; filterData(thisTrial, 0, kernel, 1, 0)];
            currentDiff(end,:) = nan(size(currentDiff(1,:)));
        end
        
        if isnan(currentDiff(calcIndex, 1))
            break;
        end
        
        thisDiffDistances = pdist2(currentDiff(currentPeaks,:), currentDiff, 'cosine');

        totalValues = 1 - (1 - thisValues ./ max(thisValues,[],2)) .* (1 - thisDiffDistances ./ max(thisDiffDistances,[],2));
    
        totalValues = totalValues(1,:);
        
        [peaks, peakIDs] = findpeaks(-totalValues);   
        
        if isempty(peakIDs)
            totalValues(isnan(totalValues)) = nanmax(totalValues);
        end
        
        allValues = totalValues;
   
        [allPeaks, allPeakIDs] = findpeaks(-totalValues);
        allPeaks = -allPeaks;
        
        [sortedPeaks, sortIDs] = sort(allPeaks);

        sortedPeakTimes = allPeakIDs(sortIDs);
        
        j = 1;
        while j < length(sortedPeakTimes)
            tempPeakTimes = sortedPeakTimes(j+1:end);
            repeatIndices = find(tempPeakTimes > sortedPeakTimes(j) - minReturnTime & tempPeakTimes < sortedPeakTimes(j) + minReturnTime);

            sortedPeakTimes(repeatIndices + j) = [];

            j = j + 1;
        end

%        

        sortedDistances = totalValues(sortedPeakTimes);
        
        

        if currentPeakCount > length(sortedDistances)            
            sigma = sortedDistances(length(sortedDistances));
            
        else
            sigma = sortedDistances(currentPeakCount);
            
%             
        end
        
        
        
        
        
        gaussianDistances = exp(-allValues.^2/(2*sigma^2));
      
        error = 1 - gaussianDistances;
        inCluster = gaussianDistances(1,:) >= exp(-1);
        outCluster = gaussianDistances(1,:) < exp(-1);

        smallestValue = min(min(gaussianDistances));
        if sum(outCluster(1,:) > 0) > 2
            totalError = mean(gaussianDistances(:, inCluster)) / (mean(gaussianDistances(:, outCluster)) + exp(-1));
        else
            totalError = 0;
        end


        totalError = totalError;

        peakSortedDistances = totalValues(1,sortedPeakTimes);
        peakDistances = exp(-peakSortedDistances.^2/(2*sigma^2));
        newCount = sum(peakDistances >= exp(-1));

        bestValues{end+1} = gaussianDistances(1,:);
        peakIndices{end+1} = allPeakIDs;
        peakValues{end+1} = totalValues(allPeakIDs);
        maxDistances(end+1) = totalError;
        worstDistances(end+1) = (mean(mean(gaussianDistances(:, outCluster))) + smallestValue);
        bestSigmas(end+1) = sigma;
        bestProjections(:,end+1) = currentSTD;
        bestDistances(:,end+1) = exp(-totalValues(1,:).^2/(2*sigma^2));
        neighborhoodCounts(end+1) = newCount;

        if verbose
            figure(3+i)
            clf;
            plot(sort(gaussianDistances, 'descend'));
            title(['k = ' num2str(currentPeakCount)]);
        end

        
        if currentPeakCount >= length(sortedPeakTimes)
            break;
        end
    end

    finalScore = maxDistances;
    distances = finalScore(MIN_NEIGHBOR_COUNT+1:end);
    distances(1) = 0;
    distances(end) = 0;
    [~, bestIndex] = max(distances);
    
    if verbose && ~isempty(gaussianDistances)
        figure(3+bestIndex)
        clf;
        plot(sort(gaussianDistances, 'descend'));
        title(['Best k = ' num2str(neighborhoodCounts(bestIndex))]);
    end

    if verbose && ~isempty(gaussianDistances)
        figure(3);
        clf;
        hold on;
        plot(finalScore(MIN_NEIGHBOR_COUNT+1:end) / max(finalScore(MIN_NEIGHBOR_COUNT+1:end)));
    end

    if isempty(bestDistances)
        distanceValues = nan(size(totalValues));
        localDistances = nan(size(totalValues));
        
        bestLocalProjections = nan(size(currentSTD))';
        bestSigmaValues = nan;
        bestDiffSigmaValues = nan;
        bestNeighborCounts = nan;
        distances = nan(size(totalValues));
    else
        distanceValues = bestDistances(:,bestIndex+MIN_NEIGHBOR_COUNT);
        plotIndices = find(distanceValues > exp(-1));

        if size(currentStream,2) > 3
            [~, plotBasis] = pca(currentStream', 'NumComponents', 3);
        else
            plotBasis = eye([size(currentStream,2) 3]);
        end

        if verbose
            figure(2);
            clf;
            hold on;
            colors = lines(8);
            h = plot3(currentStream*plotBasis(:,1),currentStream*plotBasis(:,2),currentStream*plotBasis(:,3));
            h.Color(4) = 0.2;
            scatter3(currentStream(plotIndices,:)*plotBasis(:,1),currentStream(plotIndices,:)*plotBasis(:,2),currentStream(plotIndices,:)*plotBasis(:,3), 32, distanceValues(plotIndices))
            xlabel('DPC1');
            ylabel('DPC2');
            zlabel('DPC3');
            colormap(jet(256));
        end

        test = 1;

        bestLocalProjections = bestProjections(:,bestIndex+MIN_NEIGHBOR_COUNT);
        bestSigmaValues = bestSigmas(bestIndex+MIN_NEIGHBOR_COUNT);
        localDistances = bestValues{bestIndex+MIN_NEIGHBOR_COUNT};
        bestDiffSigmaValues = 0;
        bestNeighborCounts = neighborhoodCounts(bestIndex+MIN_NEIGHBOR_COUNT);
        distances = gaussianDistances / worstDistances(bestIndex+MIN_NEIGHBOR_COUNT);

        end
end