function findLocalDistance()
    thisValues = pdist2(currentStream(calcIndex,:), currentStream);

    smoothSigma = 1;
    x = -ceil(3*smoothSigma):ceil(3*smoothSigma);
    kernel = -x.*exp(-x.^2/(2*smoothSigma^2))/(smoothSigma^3*sqrt(2*pi));

    currentDiff = filterData(currentStream, 0, kernel, 1, 0);
    thisDiffDistances = pdist2(currentDiff(calcIndex,:), currentDiff, 'cosine');
    sortedDiffDistances = sort(thisDiffDistances);
    diffCurvature = filterData(filterData(sortedDiffDistances, 0, kernel, 1, 0), 50);
    [~, bestDiffSigmaIndex] = min(diffCurvature);
    bestDiffSigma = sortedDiffDistances(ceil(bestDiffSigmaIndex/2));
    guassianDiff = exp(-thisDiffDistances.^2/(2*bestDiffSigma^2));
    guassianDiff(guassianDiff < 0.1) = 0.1;

    thisValues = thisValues ./ guassianDiff;

    [peaks, peakIDs] = findpeaks(-thisValues);
    peaks = -peaks;

    [sortedPeaks, sortIDs] = sort(peaks);

    sortedPeakTimes = peakIDs(sortIDs);

    j = 1;
    while j < length(sortedPeakTimes)
        tempPeakTimes = sortedPeakTimes(j+1:end);
        repeatIndices = find(tempPeakTimes > sortedPeakTimes(j) - MIN_RETURN_TIME & tempPeakTimes < sortedPeakTimes(j) + MIN_RETURN_TIME);

        sortedPeakTimes(repeatIndices + j) = [];


        j = j + 1;
    end

    if currentPeakCount >= length(sortedPeakTimes)
        break;
    end

    sortedDistances = thisValues(sortedPeakTimes);
    if i == 1
        sigma = thisValues(sortedPeakTimes(2));
    else
        sigma = thisValues(sortedPeakTimes(ceil(currentPeakCount/2+1)));
    end
    gaussianDistances = exp(-thisValues.^2/(2*sigma^2));
end