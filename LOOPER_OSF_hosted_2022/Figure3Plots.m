

% First load RNN or Monkey LOOPER data

romoAnalysisPipeline_forFig3
forValidation = 1;
monkeySetupForFig3
load('monkeyDataForFigure3')

%%

reconstructionPCABasis = pcaBasis';


% reconstructionPlotIndices = [1:trialLength trialLength*10+(1:trialLength) trialLength*20+(1:trialLength) trialLength*30+(1:trialLength) trialLength*40+(1:trialLength) trialLength*50+(1:trialLength)];

plotReconstruction;

%% Run once with SHOULD_VALDIATE = 0 and once with SHOULD_VALIDATE = 1

SHOULD_VALIDATE = 0;

clear tempIDs;
% tempIDs = 1:8;

Figure3PlotLoops

bestLoopIDs = [];
for i = 1:size(trialLoopIDs,1)
    for j = 1:6
        thisTrials = find(trialConditionIDs == j);
        
        bestLoopIDs(i,j) = mode(trialLoopIDs(i,thisTrials));
    end
end

percentCorrect = [];
for i = 1:size(trialLoopIDs,1)
    percentCorrect(i) = sum(trialLoopIDs(i,:) == bestLoopIDs(i, trialConditionIDs)) / length(trialConditionIDs);
end

overall = mean(percentCorrect)
F1 = mean(percentCorrect(10:33))
F1remember = min(percentCorrect(10:33))
F2 = mean(percentCorrect(40:60))

%% Run once with SHOULD_VALDIATE = 0 and once with SHOULD_VALIDATE = 1
SHOULD_VALIDATE = 1;

Figure3PlotLoops

validationPercentCorrect = [];
for i = 1:size(trialLoopIDs,1)
    validationPercentCorrect(i) = sum(trialLoopIDs(i,:) == bestLoopIDs(i, trialConditionIDs)) / length(trialConditionIDs);
end

overall = mean(validationPercentCorrect)
F1 = mean(validationPercentCorrect(15:36))
F1remember = min(validationPercentCorrect(15:36))
F2 = mean(validationPercentCorrect(40:60))

overallSim = [];
overallF1 = [];
overallF1remember = [];
overallF2 = [];
for j = 1:10000
    simulationIDs = randi(max(trialConditionIDs), size(trialLoopIDs));
    simulationPercentCorrect = [];
    for i = 1:size(trialLoopIDs,1)
        simulationPercentCorrect(i) = sum(simulationIDs(i,:) == bestLoopIDs(i, trialConditionIDs)) / length(trialConditionIDs);
    end

    overallSim(j) = mean(simulationPercentCorrect);
    overallF1(j) = mean(simulationPercentCorrect(10:33));
    overallF1remember(j) = min(simulationPercentCorrect(10:33));
    overallF2(j) = mean(simulationPercentCorrect(40:60));
end

overallP = sum(overallSim > overall) / length(overallSim)
overallF1 = sum(overallF1 > F1) / length(overallF1)
overallF1remember = sum(overallF1remember > F1remember) / length(overallF1remember)
overallF2 = sum(overallF2 > F2) / length(overallF2)




%% Reconstruction

decimateAmount = 10;

times = 1:900;
useF1s = [1 2 3];

saveData = saveData;

numTrial = max(saveData.TrialData);
trialLength = size(saveData.FinalStream,1) / numTrial;
plotTimes = 1:500;

if max(plotTimes) > trialLength
    plotTimes = plotTimes(1):trialLength;
end



allTrials = firingRatesAverage(:,useF1s,:,times);

reconstructTrials = [];
for i = 1:size(allTrials,2)
    for j = 1:size(allTrials,3)
        thisData = [];
        
        for k = 1:size(allTrials,1)
            thisData(k,:) = decimate(squeeze(allTrials(k,i,j,:)),decimateAmount);
        end
        
        reconstructTrials(:,i,j,:) = thisData;
    end
end
reconstructTrials = reconstructTrials(:,:,:,end-1-trialLength:end-2);


matSize = size(allRawFiringRates);
matData = permute(allRawFiringRates, [1, 4, 2, 3, 5]);
allTrials = reshape(matData, [matSize(1), matSize(4), matSize(2) * matSize(3) * matSize(5)]);

inputData = permute(allInputs, [1, 4, 2, 3, 5]);
allInputs = reshape(inputData, [1, matSize(4), matSize(2) * matSize(3) * matSize(5)]);

rawTrials = [];
for i = 1:size(allTrials,1)
    for j = 1:size(allTrials,3)
        rawTrials(i,:,j) = decimate(squeeze(allTrials(i,:,j)),decimateAmount);
    end
end
rawTrials = rawTrials(:,end-1-trialLength:end-2,:);
rawStream = reshape(rawTrials, [size(rawTrials,1), size(rawTrials,2)*size(rawTrials,3)])';


trialIndicies = repmat(plotTimes, [1, numTrial]);
trialIndicies = trialIndicies + kron((0:numTrial-1)*trialLength, ones(size(plotTimes)));

loopConditions = [];
for i = 1:saveData.BestLoopCount    
    allIDs = find(saveData.BestStateMap(:,1) == i);
    trialIDs = floor((allIDs-1)/trialLength) + 1;
    
    loopConditions(i) = mode(mod(trialIDs-1,6)+1);
end

clusterMeans = [];
clusterSTD = [];
clusterLengths = [];
meanTimes = [];
stateLookup = [];
for i = 1:size(saveData.BestLoopAssignments,1)
    loopID = saveData.BestLoopAssignments(i,1);
    phaseID = saveData.BestLoopAssignments(i,2);
    
    thisIndices = find(saveData.BestStateMap(:,1) == loopID & saveData.BestStateMap(:,2) == phaseID);
    
    trialTimes = mod(thisIndices, trialLength)+1;
    meanTimes(i,:) = mode(trialTimes);
    clusterLengths(i,:) = length(thisIndices);

    clusterMeans(i,:) = nanmean(rawStream(thisIndices, :), 1);
    clusterSTD(i,:) = nanstd(rawStream(thisIndices, :), [], 1);
    
    stateLookup(loopID, phaseID) = i;
end

reconstructedStream = [];
conditionCounts = zeros(size(reconstructTrials,2), size(reconstructTrials,3));
for i = 1:size(rawTrials,3)
    conditionID = mod((i-1),6)+1;
    
    F1ID = mod(conditionID - 1,3)+1;
    F2ID = floor((conditionID-1)/3)+1;
    
    thisIDs = (i-1)*trialLength + (1:trialLength);
    thisLoopIDs = saveData.BestStateMap(thisIDs,1);
    thisPhaseIDs = saveData.BestStateMap(thisIDs,2);
    
    thisStateIDs = [];
    for j = 1:length(thisLoopIDs)
        thisStateIDs(j) = stateLookup(thisLoopIDs(j), thisPhaseIDs(j));
    end
    
    conditionCounts(F1ID, F2ID) = conditionCounts(F1ID, F2ID) + 1;
    
    reconstructedStream(:, F1ID, F2ID, :, conditionCounts(F1ID, F2ID)) = clusterMeans(thisStateIDs,:)';
end

%% 
% useTimes = startIndex+1:size(rawConditions,4);

% trialRsquareds = [];
% for i = 1:size(reconstructedStream,2)
%     for j = 1:size(reconstructedStream,3)
%         for trialNum = 1:size(reconstructedStream,5)
%             targetData = squeeze(rawConditions(i,j,:,:,trialNum));
%             
%             for repeatCount = 1:size(reconstructedStream,6)
%                 thisReconstruction = squeeze(reconstructedStream(:,i,j,:,trialNum, repeatCount));
%                 
% 
%                 trialRsquareds(i,j,trialNum,repeatCount) = corr(thisReconstruction(:), targetData(:));
%             end
%         end
%     end
% end
% 
% meanTrialRsquared = nanmean(trialRsquareds(:))
% 
% 
% trialRsquareds = [];
% for i = 1:size(reconstructedStream,2)
%     for j = 1:size(reconstructedStream,3)
%         for trialNum = 1:size(reconstructedStream,5)
%             targetData = squeeze(reconstructTrials(:,i,j,:));
%             
%             for repeatCount = 1:size(reconstructedStream,6)
%                 thisReconstruction = squeeze(reconstructedStream(:,i,j,:,trialNum, repeatCount));
%                 
% 
%                 trialRsquareds(i,j,trialNum,repeatCount) = corr(thisReconstruction(:), targetData(:));
%             end
%         end
%     end
% end
% 
% meanAverageRsquared = nanmean(trialRsquareds(:))
% stdAverageRsquared = nanstd(trialRsquareds(:))



%% Simulate traces

originalTime = time;
originalTime = decimate(originalTime(times),decimateAmount);
originalTime = originalTime(end-1-trialLength:end-2);

originalLength = size(saveData.RawData,2)/numTrial;
trialLength = size(saveData.FinalStream,1) / numTrial;

% plotTime = originalTime((1:decimateAmount:trialLength*decimateAmount) + (originalLength - trialLength)*5);
startTime = 4.0;
[~, startIndex] = min(abs(originalTime - startTime));

allClusters = [];
for i = 1:length(trialData)
    thisStart = trialData{i}(:,startIndex);
    
    clusterDistances = [];
    finalStream = saveData.FinalStream;
    for j = 1:size(saveData.BestEmission,1)
        thisLoopPosition = saveData.BestLoopAssignments(j,:);
        thisIndices = find(ismember(saveData.BestStateMap, thisLoopPosition, 'rows'));
        
        stds = std(finalStream(thisIndices,:), [], 1);
        
        thisEmission = mean(finalStream(thisIndices,:), 1);
        
        clusterStream = (thisEmission - thisStart') ./ stds;           
        
        clusterDistances(j,:) = sum(clusterStream.^2,2);
    end
    
    [~, bestClusters] = min(clusterDistances, [], 1);

    allClusters(i) = bestClusters;    
end

simulateTime = size(rawTrials,2) - startIndex;
bestModel = saveData.BestModel;

rawConditions = [];
reconstructedStream = [];
conditionCounts = zeros(size(reconstructTrials,2), size(reconstructTrials,3));
for i = 1:size(rawTrials,3)
    conditionID = mod((i-1),6)+1;
    
    F1ID = mod(conditionID - 1,3)+1;
    F2ID = floor((conditionID-1)/3)+1;
    
%     thisIDs = (i-1)*trialLength + (1:trialLength);
%     thisLoopIDs = saveData.BestStateMap(thisIDs,1);
%     thisPhaseIDs = saveData.BestStateMap(thisIDs,2);
    
%     thisStateIDs = [];
%     for j = 1:length(thisLoopIDs)
%         thisStateIDs(j) = stateLookup(thisLoopIDs(j), thisPhaseIDs(j));
%     end
    
    conditionCounts(F1ID, F2ID) = conditionCounts(F1ID, F2ID) + 1;
    rawConditions(F1ID, F2ID,:,:,conditionCounts(F1ID, F2ID)) = rawTrials(:,:,i);
        
    for j = 1:100
        currentCluster = allClusters(i);

        for t = 2:simulateTime
            currentCluster(t) = randsample(1:size(bestModel,1), 1, true, bestModel(currentCluster(t-1),:));
        end

        reconstructedStream(:, F1ID, F2ID, :, conditionCounts(F1ID, F2ID), j) = clusterMeans(currentCluster,:)';
    end
end

%% Calc correlation

% Rsquareds = [];
% for i = 1:size(reconstructedStream,2)
%     for j = 1:size(reconstructedStream,3)
%         for k = 1:size(reconstructedStream,1)
%             for trialNum = 1:size(reconstructedStream,5)
%                 thisReconstruction = squeeze(reconstructedStream(k,i,j,:,trialNum));
%                 targetData = squeeze(reconstructTrials(k,i,j,:));
%                 
%                 Rsquareds(k,i,j,trialNum) = corr(thisReconstruction, targetData);
%             end
%         end
%     end
% end
% 
% meanRsquared = nanmean(Rsquareds(:))

SIM_TIME = 19;
useTimes = startIndex:startIndex+SIM_TIME-1;

% neuronRsquareds = [];
% for i = 1:size(reconstructedStream,2)
%     for j = 1:size(reconstructedStream,3)
%         for trialNum = 1:size(reconstructedStream,5)
%             targetData = squeeze(rawConditions(i,j,:,useTimes,trialNum));
%             
%             for repeatCount = 1:size(reconstructedStream,6)
%                 thisReconstruction = squeeze(reconstructedStream(:,i,j,1:SIM_TIME,trialNum, repeatCount));
%                 
%                 for k = 1:size(thisReconstruction,1)
%                     neuronRsquareds(k,i,j,trialNum,repeatCount) = corr(thisReconstruction(k,:)', targetData(k,:)');
%                 end
%             end
%         end
%     end
% end
% 
% meanNeuronRsquared = nanmean(neuronRsquareds(:))

trialRsquareds = [];
for i = 1:size(reconstructedStream,2)
    for j = 1:size(reconstructedStream,3)
        for trialNum = 1:size(reconstructedStream,5)
            targetData = squeeze(rawConditions(i,j,:,useTimes,trialNum));
            
            for repeatCount = 1:size(reconstructedStream,6)
                thisReconstruction = squeeze(reconstructedStream(:,i,j,1:SIM_TIME,trialNum, repeatCount));
                

                trialRsquareds(i,j,trialNum,repeatCount) = corr(thisReconstruction(:), targetData(:));
            end
        end
    end
end

medianTrialRsquared = nanmedian(trialRsquareds(:))
meanTrialRsquared = nanmean(trialRsquareds(:))
stdTrialRsquared = nanstd(trialRsquareds(:))

figureHandle = figure(9);
figureHandle.Renderer='Painters';
clf;
ksdensity(trialRsquareds(:));

trialRsquareds = [];
for i = 1:size(reconstructedStream,2)
    for j = 1:size(reconstructedStream,3)
        for trialNum = 1:size(reconstructedStream,5)
            targetData = squeeze(reconstructTrials(:,i,j,useTimes));
            
            for repeatCount = 1:size(reconstructedStream,6)
                thisReconstruction = squeeze(reconstructedStream(:,i,j,1:SIM_TIME,trialNum, repeatCount));
                

                trialRsquareds(i,j,trialNum,repeatCount) = corr(thisReconstruction(:), targetData(:));
            end
        end
    end
end

meanAverageRsquared = nanmean(trialRsquareds(:))
stdAverageRsquared = nanstd(trialRsquareds(:))

F1Start = 0;
F2Start = 3.5;
trialEnd = 5.5;
plotF1 = 1;
plotF2 = 1;
plotTrial = 3;
trialEnd = 5.5;

startTimes = 1:startIndex;

% plotNeurons = topFiring([6,15,21, 1]);
% thisStart = squeeze(reconstructTrials(:,plotF1,plotF2,startTimes));
% thisTarget = squeeze(reconstructTrials(:,plotF1,plotF2,useTimes));
% thisSimulated = squeeze(reconstructedStream(:,plotF1,plotF2,:,plotTrial, :));


% plotNeurons = topFiring([5,6,7,8]);
plotNeurons = [175, 23, 69, 31];
thisStart = squeeze(rawConditions(plotF1,plotF2,:,startTimes,1));
thisTarget = squeeze(rawConditions(plotF1,plotF2,:,useTimes,1));
thisSimulated = squeeze(reconstructedStream(:,plotF1,plotF2,:,plotTrial, :));
thisSimulated = thisSimulated(:,:,1);

totalFiringRates = abs(pcaBasis(:,1));
[~, topFiring] = sort(totalFiringRates, 'descend');



originalTime = time;
originalTime = decimate(originalTime(times),decimateAmount);
originalTime = originalTime(end-1-trialLength:end-2);

colors = lines(100);

figureHandle = figure(8);
figureHandle.Renderer='Painters';
clf;
hold on;
for i = 1:length(plotNeurons)
    plot(originalTime([startTimes useTimes(1)]), [thisStart(plotNeurons(i),:) thisTarget(plotNeurons(i),1)], 'Color', 'k');
end

for i = 1:length(plotNeurons)
    plot(originalTime(useTimes), thisTarget(plotNeurons(i),:), ':', 'Color', colors(i,:), 'LineWidth', 3);
end

for i = 1:length(plotNeurons)
    meanSimulated = mean(thisSimulated(:,1:SIM_TIME,:),3);
    stdSimulated = std(thisSimulated(:,1:SIM_TIME,:),[],3);
    
    plot(originalTime(useTimes), meanSimulated(plotNeurons(i),:), 'Color', colors(i,:));
    plotShadedCI(meanSimulated(plotNeurons(i),:), stdSimulated(plotNeurons(i),:), originalTime(useTimes), colors(i,:));
end
xlim([F2Start - 0.5, trialEnd])

figureHandle = figure(10);
figureHandle.Renderer='Painters';
clf;
ksdensity(trialRsquareds(:));

figureHandle = figure(11);
figureHandle.Renderer='Painters';
clf;

PLOT_CHANNELS = 1:4;
% plotChannels = [];

for i = 1:length(plotNeurons)
    fullTrace = [thisStart thisTarget];
    meanSimulated = mean(thisSimulated(:,1:SIM_TIME,:),3);
    stdSimulated = std(thisSimulated(:,1:SIM_TIME,:),[],3);
    subplot(length(plotNeurons),1,i);
    hold on;
    plot(originalTime([startTimes useTimes]), fullTrace(plotNeurons(i),:), 'k', 'LineWidth', 1);
    plotShadedCI(meanSimulated(plotNeurons(i),:), stdSimulated(plotNeurons(i),:), originalTime(useTimes), 'r');
    title(['Channel ' num2str(PLOT_CHANNELS(i))]);
    xlim([F2Start - 0.5, trialEnd])
    
%     plotChannels(i) = plotChannel;
end

%%

[correlation, index] = sort(trialRsquareds(:), 'descend');
[useF1, useF2, trialNum, repeatCount] = ind2sub(size(trialRsquareds), index(1));

% repeatCount = 1;
% useF1 = 2;
% useF2 = 1;
% trialNum = 3;

targetData = squeeze(rawConditions(useF1,useF2,:,useTimes,trialNum));
thisReconstruction = squeeze(reconstructedStream(:,useF1,useF2,1:SIM_TIME,trialNum, repeatCount));

clf
subplot(2,1,2)
imagesc(targetData)
subplot(1,2,1)
imagesc(targetData)
subplot(1,2,2)
imagesc(thisReconstruction)


% figureHandle = figure(10);
% figureHandle.Renderer='Painters';
% clf;
% ksdensity(neuronRsquareds(:));

% 
% 
% imagesc(originalTime,1:length(plotNeurons),squeeze(reconstructedStream(plotNeurons,plotF1,plotF2,:,plotTrial)));
% colormap(parula);
% xlim([min(originalTime), trialEnd])
% ylim([1 length(plotNeurons)]);
% yLimit = ylim;
% plot([F1Start F1Start], [yLimit(1) yLimit(2)], 'Color', 'k');
% plot([F1Start F1Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');
% plot([F2Start F2Start], [yLimit(1) yLimit(2)], 'Color', 'k');
% plot([F2Start F2Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');
% title('Reconstructed Data');
% 
% figureHandle = figure(9);
% figureHandle.Renderer='Painters';
% clf;
% hold on;
% imagesc(originalTime,1:length(plotNeurons),squeeze(reconstructTrials(plotNeurons,plotF1,plotF2,:)));
% colormap(parula);
% xlim([min(originalTime), trialEnd])
% ylim([1 length(plotNeurons)]);
% yLimit = ylim;
% plot([F1Start F1Start], [yLimit(1) yLimit(2)], 'Color', 'k');
% plot([F1Start F1Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');
% plot([F2Start F2Start], [yLimit(1) yLimit(2)], 'Color', 'k');
% plot([F2Start F2Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');
% title('Average Data');



% %% Sanity check
% 
% figureHandle = figure(2);
% figureHandle.Renderer='Painters';
% clf;
% hold on;
% for i = 1:numTrial
%     trialIDs = (0:trialLength-1)+(i-1)*trialLength+1;
%     
%     if IS_RNN
%         conditionID = floor((i-1)/10)+1;
%     else
%         conditionID = mod((i-1),6)+1;
%     end
%     
%     plotData = saveData.FinalStream(trialIDs,2:size(saveData.DataMean));
%     
%     h = plot(plotTime, plotData(:,1), lineStyles{conditionID});
% %     h = plot3(plotData(:,1), plotData(:,2), plotData(:,3), lineStyles{conditionID});
%     h.Color(4) = 0.2;
% end
% 
% 
% %% Setup for LOOPER
% 
% decimateAmount = 10;
% 
% saveData = saveData;
% 
% matSize = size(allBootstrappedFiringRates);
% matData = permute(allBootstrappedFiringRates, [1, 4, 2, 3, 5]);
% allTrials = reshape(matData, [matSize(1), matSize(4), matSize(2) * matSize(3) * matSize(5)]);
% 
% inputData = permute(allInputs, [1, 4, 2, 3, 5]);
% allInputs = reshape(inputData, [1, matSize(4), matSize(2) * matSize(3) * matSize(5)]);
% 
% 
% finalTrials = [];
% finalInputs = [];
% for i = 1:size(allTrials,1)
%     for j = 1:size(allTrials,3)
%         finalTrials(i,:,j) = decimate(squeeze(allTrials(i,:,j)),decimateAmount);
%         
%         if i == 1
%         	finalInputs(1,:,j) = decimate(squeeze(allInputs(1,:,j)),decimateAmount);
%         end
%     end
% end
% 
% rawData = convertToCell(finalTrials);
%                 
% lastEnd = 0;
% rawTrialData = [];
% for i = 1:length(rawData)
%     rawTrialData(lastEnd + (1:size(rawData{i}, 2))) = i;
% 
%     lastEnd = lastEnd + size(rawData{i}, 2);
% end
% 
% rawData = mergeData(rawData);
% 
% [tempData, trialData, procesedTrialSwitches] = preprocessData(rawData, [], 0, [], 0, rawTrialData, 0, true, saveData.PreprocessData.Smoothing, saveData.PreprocessData.ZScore, saveData.PreprocessData.DelayTime, saveData.PreprocessData.DelayCount, saveData.DataMean, saveData.DataSTD);
% badIndicies = [procesedTrialSwitches procesedTrialSwitches-1];
% badIndicies(badIndicies > size(tempData,2)) = [];
% badIndicies(badIndicies < 1) = [];
% % tempData(:,badIndicies) = [];
%     
% % Build cluster clouds
% 
% % clusterDistances = [];
% % finalStream = saveData.FinalStream;
% % for i = 1:size(saveData.BestEmission,1)
% %     thisLoopPosition = saveData.BestLoopAssignments(i,:);
% %     thisIndices = find(ismember(saveData.BestStateMap, thisLoopPosition, 'rows'));
% % 
% %     stds = std(finalStream(thisIndices,:), [], 1);
% % 
% %     thisEmission = mean(finalStream(thisIndices,:), 1);
% % 
% %     clusterStream = (thisEmission - tempData') ./ stds;           
% % 
% %     clusterDistances(i,:) = sum(clusterStream.^2,2);% ./ sqrt(length(thisIndices));
% % end
% 
% %     bestClusters = [];
% %     clusterDistances = pdist2(saveData.BestEmission, tempData');
% % 
% % [~, bestClusters] = min(clusterDistances, [], 1);
% % 
% % loopIDs = saveData.BestLoopAssignments(bestClusters,1);
% 
% clusterDistances = pdist2(tempData', saveData.FinalStream);
% 
% MIN_RETURN_TIME = 10;
% 
% 
% smoothSigma = 1;
% x = -ceil(3*smoothSigma):ceil(3*smoothSigma);
% kernel = -x.*exp(-x.^2/(2*smoothSigma^2))/(smoothSigma^3*sqrt(2*pi));
%     
% currentDiff = [];
% for j = 1:length(trialData)
%     thisTrial = trialData{j}';
%     currentDiff = [currentDiff; filterData(thisTrial, 0, kernel, 1, 0)];
%     currentDiff(end,:) = nan(size(currentDiff(1,:)));
% end
% 
% originalDiff = [];
% for j = 1:length(trialData)
%     thisTrial = saveData.TrialTimeSeries{j}';
%     originalDiff = [originalDiff; filterData(thisTrial, 0, kernel, 1, 0)];
%     originalDiff(end,:) = nan(size(originalDiff(1,:)));
% end
% 
% loopNearestIDs = [];
% loopWeights = [];
% waitHandle = parfor_progressbar(size(tempData,2), 'Finding best loop identities');
% for i = 1:size(tempData,2)
%     if ismember(i, badIndicies)
%         waitHandle.iterate(1);
%         continue;
%     end
%     
%     totalValues = clusterDistances(i,:);
% 
%     [peaks, peakIDs] = findpeaks(-totalValues);
%     peaks = -peaks;
% 
%     [sortedPeaks, sortIDs] = sort(peaks);
%     sortedPeakTimes = peakIDs(sortIDs);
%     
%     j = 1;
%     while j < length(sortedPeakTimes)
%         tempPeakTimes = sortedPeakTimes(j+1:end);
%         repeatIndices = find(tempPeakTimes > sortedPeakTimes(j) - MIN_RETURN_TIME & tempPeakTimes < sortedPeakTimes(j) + MIN_RETURN_TIME);
% 
%         sortedPeakTimes(repeatIndices + j) = [];
% 
% 
%         j = j + 1;
%     end
%     
%     currentPeakCount = 6;
%     
%     currentPeaks = peakIDs(sortIDs(1:currentPeakCount));
%     
%     stdPeaks = [currentPeaks-1 currentPeaks+1];
%     stdPeaks(stdPeaks < 1) = [];
%     stdPeaks(stdPeaks > size(tempData,1)) = [];
%     stdPeaks(ismember(stdPeaks, badIndicies)) = [];
%     stdPeaks = [currentPeaks, stdPeaks];
%     
%     currentSTD = std(tempData(:,stdPeaks)');
% 
%     currentSTD = currentSTD / sqrt(sum(currentSTD.^2));
%     currentSTD(currentSTD < 1e-4) = 1e-4;
% %     currentStream = tempData(:,i)' ./ currentSTD;
%     
%     originalStream = (saveData.TimeSeries' - tempData(:,i)') ./ currentSTD;
% 
%     thisValues = sqrt(sum(originalStream.^2,2));
%     
%     % Derivatives
%     
% %     thisDiffDistances = pdist2(currentDiff(i,:) ./ currentSTD, originalDiff ./ currentSTD, 'cosine')';
%     thisDiffDistances = sqrt(sum(((currentDiff(i,:) - originalDiff)./ currentSTD).^2,2));
%     
%     %Loop assignment
%     
%     totalValues = 1 - (1 - thisValues ./ max(thisValues)) .* (1 - thisDiffDistances ./ max(thisDiffDistances));
%     totalValues(badIndicies) = [];
%     
%     allPeaks = [];
%     allPeakIDs = [];
%     
%     [peaks, peakIDs] = findpeaks(-totalValues);       
% 
%     allPeaks = [allPeaks, peaks];
%     allPeakIDs = [allPeakIDs, peakIDs];
% 
%     allPeaks = -allPeaks;
%     
%     finalPeakCount = 10;
% 
%     [sortedPeaks, sortIDs] = sort(allPeaks);
% 
%     sortedPeakTimes = allPeakIDs(sortIDs);
%     
%     nearestPeakTimes = sortedPeakTimes(1:finalPeakCount);
%     
%     sortedDistances = totalValues(nearestPeakTimes);
%     
%     sigma = sortedDistances(end);
%     
%     loopNearestIDs(i,:) = saveData.BestStateMap(nearestPeakTimes);
%     loopWeights(i,:) = exp(-totalValues(nearestPeakTimes).^2/(2*sigma^2));
%     loopWeights(i,:) = loopWeights(i,:) ./ sum(loopWeights(i,:));
%     
% %     [~, bestClusters] = min(clusterDistances, [], 1);
% % 
% %     loopIDs = saveData.BestStateMap(bestClusters,1);
%     
%     waitHandle.iterate(1);
% end
% close(waitHandle);
% 
% loopIDs = mode(loopNearestIDs,2);
% removeIndicies = badIndicies;
% removeIndicies(removeIndicies > length(loopIDs)) = [];
% loopIDs(removeIndicies) = [];
% 
% 
% %%
% 
% % First build validation data
% 
% STATE_SMOOTH = 2;
% FLUX_CUTOFF = 3;
% MINIMUM_STATE_TIME = 5 *(2+1);
% IS_RNN = 0;
% 
% 
% 
% 
% 
% numTrial = max(saveData.TrialData);
% originalLength = size(saveData.RawData,2)/numTrial;
% trialLength = size(saveData.FinalStream,1) / numTrial;
% times = 1:900;
% 
% if max(times) > trialLength
%     times = times(1):trialLength;
% end
% 
% trialIndicies = repmat(times, [1, numTrial]);
% trialIndicies = trialIndicies + kron((0:numTrial-1)*trialLength, ones(size(times)));
% 
% % finalStream = saveData.FinalStream;
% 
% % [trialPCABasis, ~] = pca(finalStream(trialIndicies,:), 'NumComponents',3);
% % trialPCABAsis = pcaBasis;
% 
% % loopStarts = (0:numTrial-1)*trialLength+1;
% 
% % finalStream(loopStarts,:) = nan;
% 
% % loopIDs = saveData.BestLoopAssignments(bestClusters,1);
% if MINIMUM_STATE_TIME > 0
%     loopIDs = colfilt(loopIDs, [MINIMUM_STATE_TIME 1], 'sliding', @mode);
% end
% 
% if IS_RNN
%     loopOrder = [1, 6, 4, 3, 2, 5, 7, nan];
% 
%     lineStyles = {'r', 'g', 'b', 'r.', 'g.', 'b.'};
%     lineColors = {'r', 'g', 'b', 'r', 'g', 'b'};
% else
% %     loopOrder = [1, 2, 3, 4, 5, 6, nan];
%     loopOrder = [1, 6, 5, 4, 3, 2, nan];
% %     loopOrder = [4, 6, 5, 3, 1, 2, nan];
% 
%     lineStyles = {'r', 'g', 'b', 'r.', 'g.', 'b.'};
%     lineColors = {'r', 'g', 'b', 'r', 'g', 'b'};
%     
% %     decimateAmount = 5;
%     originalTime = [-2500 10000];
%     originalTime = originalTime(1)/1000 : 0.01 : originalTime(2)/1000;
% end
% 
% plotTime = originalTime((1:decimateAmount:trialLength*decimateAmount) + (originalLength - trialLength)*5);
% 
% figureHandle = figure(4);
% figureHandle.Renderer='Painters';
% clf;
% hold on;
% for i = 1:numTrial
%     trialIDs = (0:trialLength-1)+(i-1)*trialLength+1;
%     
%     if IS_RNN
%         conditionID = floor((i-1)/10)+1;
%     else
%         conditionID = mod((i-1),6)+1;
%     end
%     
%     IDs = loopIDs(trialIDs);
%     IDs(IDs == 0) = length(loopOrder);
%     
%     h = plot(plotTime, loopOrder(IDs) + normrnd(0,0.03,size(trialIDs)), lineStyles{conditionID});
%     h.Color(4) = 0.2;
%     
%     transitionCheck = colfilt(loopOrder(IDs), [1 2], 'sliding', @mean);
%     transitionCheck(transitionCheck == loopOrder(IDs)) = 0;
%     transitionPoints = find(transitionCheck);
%     transitionPoints = [transitionPoints, transitionPoints+1];
%     transitionPoints(transitionPoints > trialLength) = [];
%     
%     transitionData = nan(size(loopOrder(IDs)));
%     transitionData(transitionPoints) = loopOrder(IDs(transitionPoints));
%     
%     h = plot(plotTime, transitionData + normrnd(0,0.03,size(trialIDs)), lineColors{conditionID});
%     h.Color(4) = 0.2;
% end
% 
