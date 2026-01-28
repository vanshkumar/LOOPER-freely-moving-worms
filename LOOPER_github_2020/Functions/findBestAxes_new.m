function [bestLocalProjections, bestSigmaValues, bestDiffSigmaValues, localDistances, bestNeighborCounts, distanceRatios] = findBestAxes(dynamicsStream, calcIndex, meanDistanceMatrix, verbose)
    

    bestSigmaValues = [];
    bestNeighborCounts = [];
    bestLocalProjections = [];
    distanceRatios = [];
    bestDiffSigmaValues = [];
    localDistances = [];
    


    %% Near trajectories

    MIN_RETURN_TIME = 10;
    
    if ~exist('meanDistanceMatrix', 'var') || isempty(meanDistanceMatrix)
        meanDistanceMatrix = pdist2(dynamicsStream, dynamicsStream);
    end
    if ~exist('verbose', 'var') || isempty(verbose)
        verbose = 0;
    end

    
    test = meanDistanceMatrix(1:end,1:end);
    totalValues = test(calcIndex,:);

    [peaks, peakIDs] = findpeaks(-totalValues);
    peaks = -peaks;

    [sortedPeaks, sortIDs] = sort(peaks);

    sortedPeakTimes = peakIDs(sortIDs);

    i = 1;
    while i < length(sortedPeakTimes)
        tempPeakTimes = sortedPeakTimes(i+1:end);
        repeatIndices = find(tempPeakTimes > sortedPeakTimes(i) - MIN_RETURN_TIME & tempPeakTimes < sortedPeakTimes(i) + MIN_RETURN_TIME);

        sortedPeakTimes(repeatIndices + i) = [];


        i = i + 1;
    end

    %%

    MIN_NEIGHBOR_COUNT = 0;
    BEST_TRAJECTORY_COUNT = 1;

    currentPeakCount = 2;

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
    
    for i = 1:10
        if i == 1                
            currentPeaks = [peakIDs(sortIDs(1:currentPeakCount)), peakIDs(sortIDs(1:currentPeakCount)) + 1];
        else
            currentPeaks = peakIDs(sortIDs(1:currentPeakCount));
        end

        currentPCAWeight = 1;
        globalMean = mean(dynamicsStream);
        globalSTD = std(dynamicsStream - globalMean);
        currentMean = mean(dynamicsStream(currentPeaks,:));
        currentSTD = sum([BEST_TRAJECTORY_COUNT*globalSTD; currentPCAWeight*std(dynamicsStream(currentPeaks,:) - currentMean)]) / (currentPCAWeight + BEST_TRAJECTORY_COUNT);

        currentSTD = currentSTD / sqrt(sum(currentSTD.^2));
        currentStream = dynamicsStream ./ currentSTD;
        currentStream(isinf(currentStream)) = 0.0001;

        thisValues = pdist2(currentStream(currentPeaks,:), currentStream);

        smoothSigma = 1;
        x = -ceil(3*smoothSigma):ceil(3*smoothSigma);
        kernel = -x.*exp(-x.^2/(2*smoothSigma^2))/(smoothSigma^3*sqrt(2*pi));

        currentDiff = filterData(currentStream, 0, kernel, 1, 0);

        thisDiffDistances = pdist2(currentDiff(currentPeaks,:), currentDiff, 'cosine');


        totalValues = 1 - (1 - thisValues ./ max(thisValues,[],2)) .* (1 - thisDiffDistances ./ max(thisDiffDistances,[],2));
    
        allValues = totalValues(1,:);
        totalValues = totalValues(1,:);
        
        
        allPeaks = [];
        allPeakIDs = [];
        for j = 1
            [peaks, peakIDs] = findpeaks(-totalValues(j,:));            
            peakIDs = peakIDs + size(totalValues,2)*(j-1);
            
            allPeaks = [allPeaks, peaks];
            allPeakIDs = [allPeakIDs, peakIDs];
            
            allPeaks = -allPeaks;
        end
        
        [sortedPeaks, sortIDs] = sort(allPeaks);

        sortedPeakTimes = allPeakIDs(sortIDs);


        if currentPeakCount >= length(sortedPeakTimes)
            break;
        end
        
        sortedDistances = totalValues(1,sortedPeakTimes);

        totalErrors = [];
        sigmas = [];
        for j = 1:length(sortedPeakTimes)% min(length(sortedPeakTimes), currentPeakCount*2)            
            sigma = sortedDistances(j);
            
            gaussianDistances = exp(-allValues.^2/(2*sigma^2));
                
            error = 1 - gaussianDistances;
            inCluster = gaussianDistances(1,:) >= exp(-1);
            outCluster = gaussianDistances(1,:) < exp(-1);

            smallestValue = min(min(gaussianDistances));
            if sum(outCluster(1,:)) > 0
                totalError = median(gaussianDistances(:, inCluster)) / (mean(gaussianDistances(:, outCluster)) + exp(-1));
            else
                totalError = median(gaussianDistances(:, inCluster)) / (smallestValue + exp(-1));
            end
            
            totalErrors(j) = totalError;
            sigmas(j) = sortedDistances(j);
        end
        
        totalErrors = filterData(totalErrors, 1);
        [totalError, bestIndex] = max(totalErrors);
        
        sigma = sigmas(bestIndex);
        
        
        gaussianDistances = exp(-allValues.^2/(2*sigma^2));
                
        error = 1 - gaussianDistances;
        inCluster = gaussianDistances(1,:) >= exp(-1);
        outCluster = gaussianDistances(1,:) < exp(-1);

        smallestValue = min(min(gaussianDistances));
        if sum(outCluster(1,:)) > 0
            totalError = median(gaussianDistances(:, inCluster)) / (mean(gaussianDistances(:, outCluster)) + exp(-1));
        else
            totalError = median(gaussianDistances(:, inCluster)) / (smallestValue + exp(-1));
        end
        
        if i == 7
            test = 1;
        end
        
        diffDistances = diff(sortedDistances);


        totalError = totalError;

        peakSortedDistances = totalValues(1,sortedPeakTimes);
        peakDistances = exp(-peakSortedDistances.^2/(2*sigma^2));
        newCount = sum(peakDistances >= exp(-1));

        bestValues{end+1} = gaussianDistances(1,:);
        peakIndices{end+1} = allPeakIDs;
        peakValues{end+1} = totalValues(allPeakIDs);
        maxDistances(end+1) = totalError;

        bestSigmas(end+1) = sigma;
        bestProjections(:,end+1) = currentSTD;
        bestDistances(:,end+1) = exp(-totalValues(1,:).^2/(2*sigma^2));
        neighborhoodCounts(end+1) = currentPeakCount;


        if verbose
            figure(3+i)
            clf;
            plot(sort(gaussianDistances, 'descend'));
            title(['k = ' num2str(currentPeakCount)]);
        end

        currentPeakCount = ceil(max([currentPeakCount+1, currentPeakCount*1.2, newCount]));
        
        if currentPeakCount >= length(sortedPeakTimes)
            break;
        end
    end
    
    finalScore = maxDistances; 


    distances = finalScore(MIN_NEIGHBOR_COUNT+1:end);
    if 0%length(distances) > 2
        [~,peakIndices] = findpeaks(distances);
        if isempty(peakIndices)
            peakIndices = length(finalScore) - 1;
        end
        bestIndex = peakIndices(1);
    else
        [~, bestIndex] = max(distances);
    end
    [~, overrideIndex] = max(distances);
    if overrideIndex == 1
        bestIndex = 1;
    end
    
    if verbose
        figure(3+bestIndex)
        title(['Best k = ' num2str(neighborhoodCounts(bestIndex))]);
    end

    if verbose
        figure(3);
        clf;
        hold on;
        plot(finalScore(MIN_NEIGHBOR_COUNT+1:end) / max(finalScore(MIN_NEIGHBOR_COUNT+1:end)));

    end
    
    distanceValues = bestDistances(:,bestIndex+MIN_NEIGHBOR_COUNT);
    plotIndices = find(distanceValues > exp(-1));

    [plotBasis, ~] = pca(currentStream, 'NumComponents', min(4, size(currentStream,2)));

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
    distanceRatios = 0;
    
end