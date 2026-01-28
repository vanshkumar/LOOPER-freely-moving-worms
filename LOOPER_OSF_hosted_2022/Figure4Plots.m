
%% Use this to preprocess RNN data for LOOPER
% load('F:\Dropbox\ConservationOfAgentDynamics\WorkingMemoryTask\Experiments\workingMemoryCustomBroken\lstm\1\networkTrace.mat');

USE_VALIDATION = 1;

try
    TAIL_LENGTH = 1;

    percentCorrect = [];
    for i = 1:size(outputs,1)
        thisTargets = squeeze(targets(1,:,i));
        thisOutputs = double(squeeze(outputs(i,:)));

        DEBUG = 0;
        if DEBUG
            figure(1);
            clf;
            hold on;
            plot(thisTargets)
            plot(thisOutputs)
        end

        outputDiff = (thisTargets - thisOutputs);

        targetIndices = find(thisTargets ~= 0);
        checkIndicies = min(targetIndices) - TAIL_LENGTH:max(targetIndices) + TAIL_LENGTH;

        percentCorrect(i) = sum(outputDiff(checkIndicies) == 0) / length(checkIndicies);
    end

    missTrials = percentCorrect < 0.3;
    thisClass = classes == 0;
    sum(missTrials) / length(missTrials)
    sum(missTrials(thisClass)) / length(missTrials(thisClass))

    % From RNN

    NUM_TRIALS = 10;

    times = 1:220;

    useClasses = [0,1,2,3,4,5];

    trialData = permute(dynamics, [3, 1, 2]);

    allData = reshape(trialData, size(trialData,1),[]);

    [pcaBasis,~,~,~,explained] = pca(allData', 'NumComponents', 20);
    percentExplained = cumsum(explained);
    percentExplained = percentExplained(20)

%     pcaBasis = eye(size(trialData,1));

%     figure(1);
%     clf;
%     colors = lines(length(useClasses));
%     hold on;

    trialCounts = [];
    finalTrials = [];
    trialID = 1;
    for i = 1:length(useClasses)
        thisIndices = find(classes == useClasses(i) & missTrials == 0);

%         F1ID = mod(useClasses(i), 3) + 1;
%         F2ID = floor(useClasses(i) / 3) + 1;


        for j = (1:NUM_TRIALS)+USE_VALIDATION*NUM_TRIALS
            thisData = (squeeze(trialData(:,times,thisIndices(j)))' * pcaBasis)';
            for k = 1:size(thisData,1)
                finalTrials(k,:, trialID) = decimate(thisData(k,:),2);
            end
            
            trialID = trialID + 1;
        end

%         plot(mean(finalInputs(1,:,end-9:end),3), 'Color', colors(i,:));
    end
catch
end

%% Reconstruction

decimateAmount = 2;

times = 1:900;

app.SavedData = saveData;

numTrial = max(app.SavedData.TrialData);
trialLength = size(app.SavedData.FinalStream,1) / numTrial;
plotTimes = 1:500;

if max(plotTimes) > trialLength
    plotTimes = plotTimes(1):trialLength;
end



allTrials = finalTrials;

rawTrials = permute(allTrials, [1 2 3]);
rawTrials = rawTrials(:,end-1-trialLength:end-2,:);

rawStream = reshape(rawTrials, [size(rawTrials,1), size(rawTrials,2)*size(rawTrials,3)])';

rawTrials = permute(allTrials, [1 3 2]);
rawTrials = rawTrials(:,:,end-1-trialLength:end-2);

trialIndicies = repmat(plotTimes, [1, numTrial]);
trialIndicies = trialIndicies + kron((0:numTrial-1)*trialLength, ones(size(plotTimes)));

loopConditions = [];
reconstructTrials = [];
for i = 1:app.SavedData.BestLoopCount    
    allIDs = find(app.SavedData.BestStateMap(:,1) == i);
    trialIDs = floor((allIDs-1)/trialLength) + 1;
    
    loopConditions(i) = mode(floor((trialIDs-1)/10)+1);
    
    reconstructTrials(:,loopConditions(i),:) = mean(rawTrials(:,(loopConditions(i)-1)*10+(1:10),:),2);
end

clusterMeans = [];
clusterSTD = [];
clusterLengths = [];
meanTimes = [];
stateLookup = [];
for i = 1:size(app.SavedData.BestLoopAssignments,1)
    loopID = app.SavedData.BestLoopAssignments(i,1);
    phaseID = app.SavedData.BestLoopAssignments(i,2);
    
    thisIndices = find(app.SavedData.BestStateMap(:,1) == loopID & app.SavedData.BestStateMap(:,2) == phaseID);
    
    trialTimes = mod(thisIndices, trialLength)+1;
    meanTimes(i,:) = mode(trialTimes);
    clusterLengths(i,:) = length(thisIndices);

    clusterMeans(i,:) = nanmean(rawStream(thisIndices, :), 1);
    clusterSTD(i,:) = nanstd(rawStream(thisIndices, :), [], 1);
    
    stateLookup(loopID, phaseID) = i;
end

useF1s = [1,2,3];

reconstructedStream = [];
conditionCounts = zeros(length(useF1s), 2);
for i = 1:size(rawTrials,2)
    conditionID = floor((i-1)/10)+1;
    
    F1ID = mod(conditionID - 1,3)+1;
    F2ID = floor((conditionID-1)/3)+1;
    
    thisIDs = (i-1)*trialLength + (1:trialLength);
    thisLoopIDs = app.SavedData.BestStateMap(thisIDs,1);
    thisPhaseIDs = app.SavedData.BestStateMap(thisIDs,2);
    
    thisStateIDs = [];
    for j = 1:length(thisLoopIDs)
        thisStateIDs(j) = stateLookup(thisLoopIDs(j), thisPhaseIDs(j));
    end
    
    conditionCounts(F1ID, F2ID) = conditionCounts(F1ID, F2ID) + 1;
    
    reconstructedStream(:, conditionID, :, conditionCounts(F1ID, F2ID)) = clusterMeans(thisStateIDs,:)';
end

Rsquareds = [];
for i = 1:size(reconstructedStream,2)
%     for j = 1:size(reconstructedStream,3)
        for k = 1:size(reconstructedStream,1)
            for trialNumber = 1:size(reconstructedStream,4)
                thisReconstruction = squeeze(reconstructedStream(k,i,:,trialNumber));
                targetData = squeeze(reconstructTrials(k,i,:));
                
                Rsquareds(k,i,trialNumber) = corr(thisReconstruction, targetData);
            end
        end
%     end
end

meanRsquared = nanmean(Rsquareds(:))

trialRsquareds = [];
for i = 1:size(reconstructedStream,2)
%     for j = 1:size(reconstructedStream,3)
        for trialNumber = 1:size(reconstructedStream,4)
            thisReconstruction = squeeze(reconstructedStream(:,i,:,trialNumber));
            targetData = squeeze(reconstructTrials(:,i,:));

            trialRsquareds(i,trialNumber) = corr(thisReconstruction(:), targetData(:));
%         end
    end
end

meanTrialRsquared = nanmean(trialRsquareds(:))
stdTrialRsquared = nanstd(trialRsquareds(:))

%% Plot traces

F1Start = 25;
F2Start = 55;

plotF1 = 1;
plotF2 = 1;
plotCondition = (plotF2-1)*3 + plotF1;
plotTrial = 4;
plotNeurons = 1:size(reconstructedStream,1);

originalTime = time;
originalTime = decimate(originalTime(times),decimateAmount);
originalTime = originalTime(end-1-trialLength:end-2);

thisData = squeeze(reconstructTrials(plotNeurons,plotCondition,:));

distances = pdist(thisData, 'correlation');
distances(isnan(distances)) = 0;
heirarchy = linkage(distances, 'average');

leafOrder = optimalleaforder(heirarchy, distances);


figureHandle = figure(8);
figureHandle.Renderer='Painters';
clf;
subplot(2,1,1)
hold on;
imagesc(1:size(reconstructedStream,3),1:length(plotNeurons),squeeze(reconstructedStream(leafOrder,plotCondition,:,plotTrial)));
colormap(parula);
xlim([1, size(reconstructTrials,3)])
ylim([1 length(plotNeurons)]);
yLimit = ylim;
plot([F1Start F1Start], [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F1Start F1Start]+10, [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F2Start F2Start], [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F2Start F2Start]+10, [yLimit(1) yLimit(2)], 'Color', 'k');
title('Reconstructed Data');

subplot(2,1,2)
hold on;
imagesc(1:size(reconstructTrials,3),1:length(plotNeurons),squeeze(reconstructTrials(leafOrder(plotNeurons),plotCondition,:)));
colormap(parula);
xlim([1, size(reconstructTrials,3)])
ylim([1 length(plotNeurons)]);
yLimit = ylim;
plot([F1Start F1Start], [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F1Start F1Start]+10, [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F2Start F2Start], [yLimit(1) yLimit(2)], 'Color', 'k');
plot([F2Start F2Start]+10, [yLimit(1) yLimit(2)], 'Color', 'k');
title('Average Data');


% plotRaster = squeeze(-spikeRasters(leafOrder,plotCondition,:,plotTrial));
% plotRaster = convn(plotRaster, [1 1 1], 'same');
% plotRaster(plotRaster < -1) = -1;
% 
% figureHandle = figure(10);
% figureHandle.Renderer='Painters';
% clf;
% hold on;
% imagesc(time,1:length(plotNeurons),plotRaster);
% colormap(gray);
% xlim([min(originalTime), trialEnd])
% ylim([1 length(plotNeurons)]);
% yLimit = ylim;
% plot([F1Start F1Start], [yLimit(1) yLimit(2)], 'Color', 'k');
% plot([F1Start F1Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');
% plot([F2Start F2Start], [yLimit(1) yLimit(2)], 'Color', 'k');
% plot([F2Start F2Start]+0.5, [yLimit(1) yLimit(2)], 'Color', 'k');

%% Simulate like monkey

app.SavedData = saveData;

rawData = convertToCell(finalTrials);

lastEnd = 0;
rawTrialData = [];
for i = 1:length(rawData)
    rawTrialData(lastEnd + (1:size(rawData{i}, 2))) = i;

    lastEnd = lastEnd + size(rawData{i}, 2);
end

rawData = mergeData(rawData);

[tempData, trialData, procesedTrialSwitches] = preprocessData(rawData, [], 0, [], 0, rawTrialData, 0, true, saveData.PreprocessData.Smoothing, saveData.PreprocessData.ZScore, saveData.PreprocessData.DelayTime, saveData.PreprocessData.DelayCount, saveData.DataMean, saveData.DataSTD);
badIndicies = [procesedTrialSwitches procesedTrialSwitches+1];
badIndicies(badIndicies > size(tempData,2)) = [];
tempData(:,badIndicies) = [];

startIndex = 40;

allClusters = [];
for i = 1:length(trialData)
    thisStart = trialData{i}(:,startIndex);
    
    clusterDistances = [];
    finalStream = app.SavedData.FinalStream;
    for j = 1:size(app.SavedData.BestEmission,1)
        thisLoopPosition = app.SavedData.BestLoopAssignments(j,:);
        thisIndices = find(ismember(app.SavedData.BestStateMap, thisLoopPosition, 'rows'));
        
        stds = std(finalStream(thisIndices,:), [], 1);
        
        thisEmission = mean(finalStream(thisIndices,:), 1);
        
        clusterStream = (thisEmission - thisStart') ./ stds;           
        
        clusterDistances(j,:) = sum(clusterStream.^2,2);
    end
    
    [~, bestClusters] = min(clusterDistances, [], 1);

    allClusters(i) = bestClusters;        
    
end

simulateTime = size(rawTrials,3) - startIndex;
bestModel = app.SavedData.BestModel;

rawConditions = [];
reconstructedStream = [];
conditionCounts = zeros(size(reconstructTrials,2), size(reconstructTrials,3));
for i = 1:length(trialData)
    conditionID = floor((i-1)/10)+1;
    
    F1ID = mod(conditionID - 1,3)+1;
    F2ID = floor((conditionID-1)/3)+1;
    
%     thisIDs = (i-1)*trialLength + (1:trialLength);
%     thisLoopIDs = app.SavedData.BestStateMap(thisIDs,1);
%     thisPhaseIDs = app.SavedData.BestStateMap(thisIDs,2);
    
%     thisStateIDs = [];
%     for j = 1:length(thisLoopIDs)
%         thisStateIDs(j) = stateLookup(thisLoopIDs(j), thisPhaseIDs(j));
%     end
    
    conditionCounts(conditionID) = conditionCounts(conditionID) + 1;
    rawConditions(conditionID,:,:,conditionCounts(conditionID)) = squeeze(rawTrials(:,i,:));
        
    for j = 1:100
        currentCluster = allClusters(i);

        for t = 2:simulateTime
            currentCluster(t) = randsample(1:size(bestModel,1), 1, true, bestModel(currentCluster(t-1),:));
        end

        reconstructedStream(:, conditionID, :, conditionCounts(conditionID), j) = clusterMeans(currentCluster,:)';
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

SIM_TIME = 40;
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
%     for j = 1:size(reconstructedStream,3)
        for trialNum = 1:size(reconstructedStream,4)
            targetData = squeeze(rawConditions(i,:,useTimes,trialNum));
            
            for repeatCount = 1:size(reconstructedStream,5)
                thisReconstruction = squeeze(reconstructedStream(:,i,1:SIM_TIME,trialNum, repeatCount));
                

                trialRsquareds(i,trialNum,repeatCount) = corr(thisReconstruction(:), targetData(:));
            end
        end
%     end
end

medianTrialRsquared = nanmedian(trialRsquareds(:))
meanTrialRsquared = nanmean(trialRsquareds(:))
stdTrialRsquared = nanstd(trialRsquareds(:))


trialRsquareds = [];
for i = 1:size(reconstructedStream,2)
%     for j = 1:size(reconstructedStream,3)
        for trialNum = 1:size(reconstructedStream,4)
            targetData = squeeze(reconstructTrials(:,i,useTimes));
            
            for repeatCount = 1:size(reconstructedStream,5)
                thisReconstruction = squeeze(reconstructedStream(:,i,1:SIM_TIME,trialNum, repeatCount));
                

                trialRsquareds(i,trialNum,repeatCount) = corr(thisReconstruction(:), targetData(:));
            end
        end
%     end
end

meanAverageRsquared = nanmean(trialRsquareds(:))
stdAverageRsquared = nanstd(trialRsquareds(:))

F1Start = 25;
F2Start = 35;
plotCondition = 1;
plotTrial = 1;
trialEnd = 88;

startTimes = 1:startIndex;

thisStart = pcaBasis * squeeze(reconstructTrials(:,plotCondition,startTimes));
thisTarget = pcaBasis * squeeze(reconstructTrials(:,plotCondition,useTimes));
thisSimulated = pcaBasis * reshape(squeeze(reconstructedStream(:,plotCondition,:,plotTrial, :)), size(reconstructTrials,1), []);

thisSimulated = reshape(thisSimulated, size(thisSimulated,1), size(reconstructedStream,3), size(reconstructedStream,5));

totalFiringRates = abs(pcaBasis(:,1));
[~, topFiring] = sort(totalFiringRates, 'descend');

plotNeurons = topFiring([1,2,31,50]);

originalTime = 1:size(reconstructTrials,3);
% originalTime = time;
% originalTime = decimate(originalTime(times),decimateAmount);
% originalTime = originalTime(end-1-trialLength:end-2);

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
xlim([plotTrial, trialEnd])

figureHandle = figure(9);
figureHandle.Renderer='Painters';
clf;
ksdensity(trialRsquareds(:));

figureHandle = figure(10);
figureHandle.Renderer='Painters';
clf;

for i = 1:length(plotNeurons)
    fullTrace = [thisStart thisTarget];
    meanSimulated = mean(thisSimulated(:,1:SIM_TIME,:),3);
    stdSimulated = std(thisSimulated(:,1:SIM_TIME,:),[],3);
    subplot(length(plotNeurons),1,i);
    hold on;
    plot(originalTime([startTimes useTimes]), fullTrace(plotNeurons(i),:), 'k', 'LineWidth', 1);
    plotShadedCI(meanSimulated(plotNeurons(i),:), stdSimulated(plotNeurons(i),:), originalTime(useTimes), 'r');
    title(['Channel ' num2str(plotChannel)]);
    xlim([F2Start - 10, trialEnd])
    
    plotChannels(i) = plotChannel;
end


%% Good RNN
load('F:\Dropbox\ConservationOfAgentDynamics\WorkingMemoryTask\Experiments\workingMemoryCustomFrequenciesLong\lstm\1\networkTrace.mat');
load('goodRNNResult');

figureNumber = 1;
IS_VALIDATION = 0;

Figure4PlotLoops

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
F1 = mean(percentCorrect(14:35))
F1remember = min(percentCorrect(14:35))
F2 = mean(percentCorrect(40:60))

%% Reconstruction

reconstructionPCABasis = pcaBasis';

reconstructionPlotIndices = [1:trialLength trialLength*10+(1:trialLength) trialLength*20+(1:trialLength) trialLength*30+(1:trialLength) trialLength*40+(1:trialLength) trialLength*50+(1:trialLength)];

plotReconstruction;



%% Validation

IS_VALIDATION = 1;
figureNumber = 2;

Figure4PlotLoops

validationPercentCorrect = [];
for i = 1:size(trialLoopIDs,1)
    validationPercentCorrect(i) = sum(trialLoopIDs(i,:) == bestLoopIDs(i, trialConditionIDs)) / length(trialConditionIDs);
end

overall = mean(percentCorrect)
F1 = mean(percentCorrect(14:35))
F1remember = min(percentCorrect(14:35))
F2 = mean(percentCorrect(40:60))


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
    overallF1(j) = mean(simulationPercentCorrect(14:35));
    overallF1remember(j) = min(simulationPercentCorrect(14:35));
    overallF2(j) = mean(simulationPercentCorrect(40:60));
end

overallP = sum(overallSim > overall) / length(overallSim)
overallF1 = sum(overallF1 > F1) / length(overallF1)
overallF1remember = sum(overallF1remember > F1remember) / length(overallF1remember)
overallF2 = sum(overallF2 > F2) / length(overallF2)


%% Broken RNN
load('F:\Dropbox\ConservationOfAgentDynamics\WorkingMemoryTask\Experiments\workingMemoryCustomBroken\lstm\1\networkTrace.mat');
load('rnnBroken1ForFig4');

IS_VALIDATION = 0;
figureNumber = 3;

Figure4PlotLoops

%% Find miss trials
PLOT_BROKEN = 1;

RNNAccuracies = {};

for PLOT_BROKEN = 0:1

    if ~PLOT_BROKEN
        % Good RNN
        load('F:\Dropbox\ConservationOfAgentDynamics\WorkingMemoryTask\Experiments\workingMemoryCustomFrequenciesLong\lstm\1\networkTraceTests.mat')
        figureNumber = 4;
    else
        % Broken RNN
        load('F:\Dropbox\ConservationOfAgentDynamics\WorkingMemoryTask\Experiments\workingMemoryCustomBroken\lstm\1\networkTraceTests.mat')
        figureNumber = 5;
    end

    TAIL_LENGTH = 0;

    allOutputs = [];
    allTargets = [];
    percentCorrect = [];
    for i = 1:size(outputs,1)
        thisTargets = squeeze(targets(1,:,i));
        thisOutputs = double(squeeze(outputs(i,:)));

        DEBUG = 0;
        if DEBUG

            figure(1);
            clf;
            hold on;
            plot(thisTargets)
            plot(thisOutputs)
        end


        outputDiff = (thisTargets - thisOutputs);

        targetIndices = find(thisTargets ~= 0);
        checkIndicies = min(targetIndices) - TAIL_LENGTH:max(targetIndices) + TAIL_LENGTH;


        allOutputs(i,:) = thisOutputs(checkIndicies);
        allTargets(i,:) = thisTargets(checkIndicies);

        percentCorrect(i) = sum(outputDiff(checkIndicies) == 0) / length(checkIndicies);
    end

    missTrials = percentCorrect < 0.3;
    hitTrials = percentCorrect > 0.5;

    sum(missTrials) / length(missTrials)

    maxClass = max(classes);

    accuracies = [];
    hitAccuracy = [];
    targetOutputs = [];
    for i = 1:maxClass+1
        thisIDs = find(classes == i-1);

    %     hitAccuracy(i) = mean(percentCorrect(thisIDs));
        thisOutputs = allOutputs(thisIDs,:);

        thisTargets = allTargets(thisIDs,:);
        targetOutputs(i) = mode(thisTargets(:));

        hitAccuracy(i,:) = histcounts(thisOutputs(:), [-0.5, 0.5, 1.5, 2.5], 'Normalization', 'pdf');

        for j = 1:2        
            if j == 1
                accuracies(i,j) = length(find(classes == i-1 & hitTrials));
            else
                accuracies(i,j) = length(find(classes == i-1 & missTrials));
            end
        end
    end

    missConditions = find(sum(accuracies >= 10,2) >= 2);

    figureHandle = figure(figureNumber);
    figureHandle.Renderer='Painters';
    clf
    hold on;
    imagesc(hitAccuracy);
    yticks(1:length(hitAccuracy));

    if length(hitAccuracy) > 15
        yticklabels({   '20Hz -> 5Hz', ...
                        '20Hz -> 15Hz', ...
                        '20Hz -> 25Hz', ...
                        '20Hz -> 30Hz', ...
                        '20Hz -> 50Hz', ...
                        '40Hz -> 5Hz', ...
                        '40Hz -> 15Hz', ...
                        '40Hz -> 25Hz', ...
                        '40Hz -> 30Hz', ...
                        '40Hz -> 50Hz', ...
                        '10Hz -> 5Hz', ...
                        '10Hz -> 15Hz', ...
                        '10Hz -> 25Hz', ...
                        '10Hz -> 30Hz', ...
                        '10Hz -> 50Hz', ...
                        '0Hz -> 5Hz', ...
                        '0Hz -> 15Hz', ...
                        '0Hz -> 25Hz', ...
                        '0Hz -> 30Hz', ...
                        '0Hz -> 50Hz'});
    else
        yticklabels({   '30Hz -> 20Hz', ...
                        '30Hz -> 30Hz', ...
                        '30Hz -> 40Hz', ...
                        '30Hz -> 50Hz', ...
                        '30Hz -> 60Hz', ...
                        '40Hz -> 20Hz', ...
                        '40Hz -> 30Hz', ...
                        '40Hz -> 40Hz', ...
                        '40Hz -> 50Hz', ...
                        '40Hz -> 60Hz', ...
                        '50Hz -> 20Hz', ...
                        '50Hz -> 30Hz', ...
                        '50Hz -> 40Hz', ...
                        '50Hz -> 50Hz', ...
                        '50Hz -> 60Hz'});
    end
    xticks(1:3);
    xticklabels({'No output', 'Less than', 'Greater than'});
    scatter(targetOutputs + 1, 1:length(targetOutputs), 'rx');
    ylim([0.5 length(targetOutputs) + 0.5]);
    
    RNNAccuracies{end+1} = hitAccuracy;
end

%% Reconstruction

load('rnnBroken1ForSimulation');

app.SavedData = saveData;

numTrial = max(app.SavedData.TrialData);
trialLength = size(app.SavedData.FinalStream,1) / numTrial;
trialStarts = 1:trialLength:size(app.SavedData.FinalStream,1);

stateModel = app.SavedData.BestModel;

stateMap = [];
for i = 1:size(stateModel,1)
    loopID = app.SavedData.BestLoopAssignments(i,1);
    phaseID = app.SavedData.BestLoopAssignments(i,2);
    
    stateIndices = find(app.SavedData.BestStateMap(:,1) == loopID & app.SavedData.BestStateMap(:,2) == phaseID);

    stateMap(stateIndices) = i;
end

% stateEmissions = [];
% stateSTDs = [];
% 
% stateEmissions(i,:) = mean(app.SavedData.FinalStream(stateIndices,1:20),1);
% stateSTDs(i,:) = std(app.SavedData.FinalStream(stateIndices,1:20),1);


rawData = app.SavedData.RawData;
rawInputs = app.SavedData.Inputs;
rawOutputs = app.SavedData.Outputs;

trialStates = [];
trialActivity = [];
trialInputs = [];
trialOutputs = [];
simulationInputs = [];
simulationTargets = [];
trueActivity = [];
for i = 1:numTrial
    thisIndices = find(app.SavedData.TrialData == i);
    finalindices = trialStarts(i):trialStarts(i)+trialLength-1;
    
    thisData = rawData(:, thisIndices);
    thisInputs = rawInputs(:, thisIndices);
    thisOutputs = rawOutputs(:, thisIndices);
    
    finalTimes = size(thisData,2)-trialLength+1-2:size(thisData,2)-2;
    
    trialActivity = [trialActivity thisData(:,finalTimes)];
%     trialInputs = [trialInputs finalInputs(:,finalTimes,i)];
%     trialOutputs = [trialOutputs finalOutputs(:,finalTimes,i)];
    trialInputs = [trialInputs thisInputs(:,finalTimes)];
    trialOutputs = [trialOutputs thisOutputs(:,finalTimes)];
    trialStates = [trialStates stateMap(finalindices)];
    
%     simulationInputs(:,:,i) = finalInputs(:,finalTimes,i);
%     simulationTargets(:,:,i) = finalTargets(:,finalTimes,i);
    simulationInputs(:,:,i) = thisInputs(:,finalTimes);
    simulationTargets(:,:,i) = thisOutputs(:,finalTimes);
    trueActivity(:,:,i) = thisData(:,finalTimes);
end

stateOuputDistributions = [];
stateInputMeans = [];
stateInputSTDs = [];
stateMeans = [];
for i = 1:size(stateModel,1)
    loopID = app.SavedData.BestLoopAssignments(i,1);
    phaseID = app.SavedData.BestLoopAssignments(i,2);
    
    stateIndices = find(app.SavedData.BestStateMap(:,1) == loopID & app.SavedData.BestStateMap(:,2) == phaseID);

    stateOuputDistributions(i,:) = histcounts(trialOutputs(stateIndices), [-0.5 0.5 1.5 2.5], 'Normalization', 'pdf');
    stateInputMeans(i) = mean(trialInputs(stateIndices));
    
    networkNoise = 2;
    stateInputSTDs(i) = max(networkNoise, std(trialInputs(stateIndices)));
    
    stateMeans(i,:) = mean(trialActivity(:,stateIndices),2);
end

startState = mode(trialStates(trialStarts));


%% Simulate model

F1Start = 25;
F1End = 29;
F2Start = 55;
F2End = 59;

% F1s = [10 40 20 10 40 20];
% F2s = [5 15 30 15 25 50];
% targets = [1 1 2 2 1 2];
F1s = [10 20 40 10 20 40];
F2s = [5 15 30 15 25 50];
targets = [1 1 1 2 2 2];

SIMULATE_COUNT = 100;

inputs = simulationInputs;
% targets = simulationTargets;

simulationBackTimes = [];
simulatedTargets = [];
simualtedTraces = [];
averageActivities = [];
for i = 1:6%size(inputs,3)
%     thisInput = filterData(inputs(:,:,i),1);
    thisTrials = (i-1)*10+(1:10);
    averageActivities(:,:,i) = mean(trueActivity(:,:,thisTrials),3);
%     thisInput = mean(inputs(:,:,thisTrials),3);
%     thisInput = colfilt(thisInput, [1 5], 'sliding', @max);
    thisInput = zeros(size(inputs(:,:,1)));
    thisInput(F1Start:F1End) = F1s(i);
    thisInput(F2Start:F2End) = F2s(i);    
%     thisInput = colfilt(thisInput, [1 5], 'sliding', @max);

    stepCounter = 1;
    
    for j = 1:SIMULATE_COUNT
        currentState(1) = startState;    
        simulatedOutputs(1) = randsample(0:2, 1, true, stateOuputDistributions(currentState(1),:));
        simualtedTraces(:,1,i,j) = stateMeans(currentState(1),:);

        for t = 2:size(inputs,2)
            totalWeight = 0;
            useMatrix = stateModel;
            backTime = 1;
            while totalWeight == 0
                stateWeights = useMatrix(currentState(t-backTime),:);
                inputWeights = exp(-(thisInput(t) - stateInputMeans).^2 ./ (2*stateInputSTDs.^2));
                inputWeights(isnan(inputWeights)) = 0;
                inputWeights(inputWeights < exp(-2)) = 0;

                stateWeights = stateWeights .* inputWeights;
                totalWeight = sum(stateWeights);
                
                if totalWeight > 0
                    break;
                end
                
                useMatrix = useMatrix * stateModel;
                backTime = backTime + 1;
            end
            stateWeights = stateWeights / sum(stateWeights);
            
            simulationBackTimes(i,stepCounter) = backTime;
            stepCounter = stepCounter + 1;

            currentState(t) = randsample(1:size(stateModel,1), 1, true, stateWeights);
            simulatedOutputs(t) = randsample(0:2, 1, true, stateOuputDistributions(currentState(t),:));
            simulatedTraces(:,t,i,j) = stateMeans(currentState(t),:);
%             disp(['t ' num2str(t)]);
%             disp(['currentState ' num2str(currentState(t))]);
%             disp(['stateInputMeans ' num2str(stateInputMeans(currentState(t)))]);
%             disp(['thisInput ' num2str(thisInput(t))]);
        end
        
        simulatedTargets(:,i,j) = simulatedOutputs;
    end
end

figureHandle = figure(6);
figureHandle.Renderer='Painters';
clf;
imagesc(squeeze(simulatedTargets(:,1,:))');
ylabel('Tiral ID')
xlabel('Time')
% figure(2);
% clf;
% imagesc(squeeze(simulatedTargets(:,2,:)));

trialAccuracy = [];
for i = 1:size(simulatedTargets,2)
    hits = [];
    for j = 1:size(simulatedTargets,3)
        trueTarget = targets(i);
        
        thisOutputs = simulatedTargets(60:65,i,j)';
        thisOutputs(thisOutputs == 0) = [];
        thisOutput = mode(thisOutputs);
        
        hits(j) = (thisOutput == trueTarget);
    end
    
    trialAccuracy(i) = sum(hits) / size(simulatedTargets,3);
end
trialAccuracy

RSquareds = [];
for i = 1:size(simulatedTraces,3)
    trueAverage = averageActivities(:,:,i);
        
    for j = 1:size(simulatedTraces,4)
        simData = simulatedTraces(:,:,i,j);
        
        RSquareds(i,j) = corr(trueAverage(:), simData(:));
    end
end
meanRSquared = mean(RSquareds,2)
totalMeanRSquared = mean(meanRSquared)

allTraces = squeeze(simulatedTraces(1,:,1,:));
distances = pdist2(squeeze(averageActivities(1,:,1)), allTraces', 'correlation');

[~, bestTrial] = min(distances);

figureHandle = figure(6);
figureHandle.Renderer='Painters';
clf;
subplot(2,1,1);
imagesc(averageActivities(:,:,1));
ylabel('PC #');
title('Average')
xlabel('Time')
subplot(2,1,2);
imagesc(simulatedTraces(:,:,1,bestTrial));
title('Simulated')
ylabel('PC #');
xlabel('Time')

%% Simualte bad trials

F1s = [20, 20, 20, 20, 20, 40, 40, 40, 40, 40, 10, 10, 10, 10, 10];
F2s = [ 5, 15, 25, 30, 50,  5, 15, 25, 30, 50,  5, 15, 25, 30, 50];
targets = [1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 1, 2, 2, 2, 2];

yticklabels({   '20Hz -> 5Hz', ...
                '20Hz -> 15Hz', ...
                '20Hz -> 25Hz', ...
                '20Hz -> 30Hz', ...
                '20Hz -> 50Hz', ...
                '40Hz -> 5Hz', ...
                '40Hz -> 15Hz', ...
                '40Hz -> 25Hz', ...
                '40Hz -> 30Hz', ...
                '40Hz -> 50Hz', ...
                '10Hz -> 5Hz', ...
                '10Hz -> 15Hz', ...
                '10Hz -> 25Hz', ...
                '10Hz -> 30Hz', ...
                '10Hz -> 50Hz'});

SIMULATE_COUNT = 100;

inputs = simulationInputs;
% targets = simulationTargets;

simulationBackTimes = [];
simulatedTargets = [];
simualtedTraces = [];
averageActivities = [];
for i = 1:size(F1s,2)%size(inputs,3)
%     thisInput = filterData(inputs(:,:,i),1);
    thisTrials = (i-1)*10+(1:10);
%     averageActivities(:,:,i) = mean(trueActivity(:,:,thisTrials),3);
%     thisInput = mean(inputs(:,:,thisTrials),3);
%     thisInput = colfilt(thisInput, [1 5], 'sliding', @max);
    thisInput = zeros(size(inputs(:,:,1)));
    thisInput(F1Start:F1End) = F1s(i);
    thisInput(F2Start:F2End) = F2s(i);    
%     thisInput = colfilt(thisInput, [1 5], 'sliding', @max);

    stepCounter = 1;
    
    for j = 1:SIMULATE_COUNT
        currentState(1) = startState;    
        simulatedOutputs(1) = randsample(0:2, 1, true, stateOuputDistributions(currentState(1),:));
        simualtedTraces(:,1,i,j) = stateMeans(currentState(1),:);

        for t = 2:size(inputs,2)
            totalWeight = 0;
            useMatrix = stateModel;
            backTime = 1;
            while totalWeight == 0
                stateWeights = useMatrix(currentState(t-backTime),:);
                inputWeights = exp(-(thisInput(t) - stateInputMeans).^2 ./ (2*stateInputSTDs.^2));
                inputWeights(isnan(inputWeights)) = 0;
                inputWeights(inputWeights < exp(-2)) = 0;

                stateWeights = stateWeights .* inputWeights;
                totalWeight = sum(stateWeights);
                
                if totalWeight > 0
                    break;
                end
                
                useMatrix = useMatrix * stateModel;
                backTime = backTime + 1;
            end
            stateWeights = stateWeights / sum(stateWeights);
            
            simulationBackTimes(i,stepCounter) = backTime;
            stepCounter = stepCounter + 1;

            currentState(t) = randsample(1:size(stateModel,1), 1, true, stateWeights);
            simulatedOutputs(t) = randsample(0:2, 1, true, stateOuputDistributions(currentState(t),:));
            simulatedTraces(:,t,i,j) = stateMeans(currentState(t),:);
%             disp(['t ' num2str(t)]);
%             disp(['currentState ' num2str(currentState(t))]);
%             disp(['stateInputMeans ' num2str(stateInputMeans(currentState(t)))]);
%             disp(['thisInput ' num2str(thisInput(t))]);
        end
        
        simulatedTargets(:,i,j) = simulatedOutputs;
    end
end

trialAccuracy = [];
for i = 1:size(simulatedTargets,2)
    hits = [];
    for j = 1:size(simulatedTargets,3)
        trueTarget = targets(i);
        
        thisOutputs = simulatedTargets(60:65,i,j)';
        thisOutputs(thisOutputs == 0) = [];
        thisOutput = mode(thisOutputs);
        
        hits(j) = (thisOutput == trueTarget);
    end
    
    trialAccuracy(i) = sum(hits) / size(simulatedTargets,3);
end
trialAccuracy

%% Display all accuracies

% 1 '30Hz -> 20Hz', ...
% 2 '30Hz -> 30Hz', ...
% 3 '30Hz -> 40Hz', ...
% 4 '30Hz -> 50Hz', ...
% 5 '30Hz -> 60Hz', ...
% 6 '40Hz -> 20Hz', ...
% 7 '40Hz -> 30Hz', ...
% 8 '40Hz -> 40Hz', ...
% 9 '40Hz -> 50Hz', ...
% 10 '40Hz -> 60Hz', ...
% 11 '50Hz -> 20Hz', ...
% 12 '50Hz -> 30Hz', ...
% 13 '50Hz -> 40Hz', ...
% 14 '50Hz -> 50Hz', ...
% 15 '50Hz -> 60Hz'});
goodRNN = RNNAccuracies{1}([15 13 3 1 13 11 5 3],:);
goodRNNHits = [2 1 2 1 1 1 2 2];
goodRNNAccuracy = [];
for i = 1:size(goodRNN,1)
    goodRNNAccuracy(i) = goodRNN(i,goodRNNHits(i)+1);
end

% 1 '20Hz -> 5Hz', ...
% 2 '20Hz -> 15Hz', ...
% 3 '20Hz -> 25Hz', ...
% 4 '20Hz -> 30Hz', ...
% 5 '20Hz -> 50Hz', ...
% 6 '40Hz -> 5Hz', ...
% 7 '40Hz -> 15Hz', ...
% 8 '40Hz -> 25Hz', ...
% 9 '40Hz -> 30Hz', ...
% 10 '40Hz -> 50Hz', ...
% 11 '10Hz -> 5Hz', ...
% 12 '10Hz -> 15Hz', ...
% 13 '10Hz -> 25Hz', ...
% 14 '10Hz -> 30Hz', ...
% 15 '10Hz -> 50Hz'
badRNN = RNNAccuracies{2}([10 9 3 2 8 7 5 4],:);
badRNNHits = [2 1 2 1 1 1 2 2];
badRNNAccuracy = [];
for i = 1:size(goodRNN,1)
    badRNNAccuracy(i) = badRNN(i,badRNNHits(i)+1);
end

simulationAccuracy = trialAccuracy([10 9 3 2 8 7 5 4]);

allAccuracies = [goodRNNAccuracy; badRNNAccuracy; simulationAccuracy];

figureHandle = figure(7);
figureHandle.Renderer='Painters';
clf;
imagesc(allAccuracies');

yticklabels({   '7 -> 8', ...
                '7 -> 6', ...
                '4 -> 5', ...
                '4 -> 3', ...
                '7 -> 5', ...
                '7 -> 3', ...
                '4 -> 8', ...
                '4 -> 6'});


%%

% %%
% 
% 
% 
% 
% 
% 
% decimateAmount = 2;
% 
% times = 1:900;
% usF1s = [1 2 3];
% 
% app.SavedData = saveData;
% 
% numTrial = max(app.SavedData.TrialData);
% trialLength = size(app.SavedData.FinalStream,1) / numTrial;
% plotTimes = 1:500;
% 
% if max(plotTimes) > trialLength
%     plotTimes = plotTimes(1):trialLength;
% end
% 
% 
% 
% allTrials = firingRatesAverage(:,useF1s,:,times);
% 
% reconstructTrials = [];
% for i = 1:size(allTrials,2)
%     for j = 1:size(allTrials,3)
%         thisData = [];
%         
%         for k = 1:size(allTrials,1)
%             thisData(k,:) = decimate(squeeze(allTrials(k,i,j,:)),decimateAmount);
%         end
%         
%         reconstructTrials(:,i,j,:) = thisData;
%     end
% end
% reconstructTrials = reconstructTrials(:,:,:,end-1-trialLength:end-2);
% 
% 
% matSize = size(allRawFiringRates);
% matData = permute(allRawFiringRates, [1, 4, 2, 3, 5]);
% allTrials = reshape(matData, [matSize(1), matSize(4), matSize(2) * matSize(3) * matSize(5)]);
% 
% inputData = permute(allInputs, [1, 4, 2, 3, 5]);
% allInputs = reshape(inputData, [1, matSize(4), matSize(2) * matSize(3) * matSize(5)]);
% 
% rawTrials = [];
% for i = 1:size(allTrials,1)
%     for j = 1:size(allTrials,3)
%         rawTrials(i,:,j) = decimate(squeeze(allTrials(i,:,j)),decimateAmount);
%     end
% end
% rawTrials = rawTrials(:,end-1-trialLength:end-2,:);
% rawStream = reshape(rawTrials, [size(rawTrials,1), size(rawTrials,2)*size(rawTrials,3)])';
% 
% 
% trialIndicies = repmat(plotTimes, [1, numTrial]);
% trialIndicies = trialIndicies + kron((0:numTrial-1)*trialLength, ones(size(plotTimes)));
% 
% loopConditions = [];
% for i = 1:app.SavedData.BestLoopCount    
%     allIDs = find(app.SavedData.BestStateMap(:,1) == i);
%     trialIDs = floor((allIDs-1)/trialLength) + 1;
%     
%     loopConditions(i) = mode(mod(trialIDs-1,6)+1);
% end
% 
% clusterMeans = [];
% clusterSTD = [];
% clusterLengths = [];
% meanTimes = [];
% stateLookup = [];
% for i = 1:size(app.SavedData.BestLoopAssignments,1)
%     loopID = app.SavedData.BestLoopAssignments(i,1);
%     phaseID = app.SavedData.BestLoopAssignments(i,2);
%     
%     thisIndices = find(app.SavedData.BestStateMap(:,1) == loopID & app.SavedData.BestStateMap(:,2) == phaseID);
%     
%     trialTimes = mod(thisIndices, trialLength)+1;
%     meanTimes(i,:) = mode(trialTimes);
%     clusterLengths(i,:) = length(thisIndices);
% 
%     clusterMeans(i,:) = nanmean(rawStream(thisIndices, :), 1);
%     clusterSTD(i,:) = nanstd(rawStream(thisIndices, :), [], 1);
%     
%     stateLookup(loopID, phaseID) = i;
% end
% 
% % reconstructedStream = [];
% % conditionCounts = zeros(size(reconstructTrials,2), size(reconstructTrials,3));
% % for i = 1:size(rawTrials,3)
% %     conditionID = mod((i-1),6)+1;
% %     
% %     F1ID = mod(conditionID - 1,3)+1;
% %     F2ID = floor((conditionID-1)/3)+1;
% %     
% %     thisIDs = (i-1)*trialLength + (1:trialLength);
% %     thisLoopIDs = app.SavedData.BestStateMap(thisIDs,1);
% %     thisPhaseIDs = app.SavedData.BestStateMap(thisIDs,2);
% %     
% %     thisStateIDs = [];
% %     for j = 1:length(thisLoopIDs)
% %         thisStateIDs(j) = stateLookup(thisLoopIDs(j), thisPhaseIDs(j));
% %     end
% %     
% %     conditionCounts(F1ID, F2ID) = conditionCounts(F1ID, F2ID) + 1;
% %     
% %     reconstructedStream(:, F1ID, F2ID, :, conditionCounts(F1ID, F2ID)) = clusterMeans(thisStateIDs,:)';
% % end
% 
% 
% %% Simulate traces
% 
% startTime = 3.8;
% [~, startIndex] = min(abs(plotTime - startTime));
% 
% allClusters = [];
% for i = 1:length(trialData)
%     thisStart = trialData{i}(:,startIndex);
%     
%     clusterDistances = [];
%     finalStream = app.SavedData.FinalStream;
%     for j = 1:size(app.SavedData.BestEmission,1)
%         thisLoopPosition = app.SavedData.BestLoopAssignments(j,:);
%         thisIndices = find(ismember(app.SavedData.BestStateMap, thisLoopPosition, 'rows'));
%         
%         stds = std(finalStream(thisIndices,:), [], 1);
%         
%         thisEmission = mean(finalStream(thisIndices,:), 1);
%         
%         clusterStream = (thisEmission - thisStart') ./ stds;           
%         
%         clusterDistances(j,:) = sum(clusterStream.^2,2);
%     end
%     
%     [~, bestClusters] = min(clusterDistances, [], 1);
% 
%     allClusters(i) = bestClusters;        
%     
% end
% 
% simulateTime = size(rawTrials,2) - startIndex;
% bestModel = app.SavedData.BestModel;
% 
% rawConditions = [];
% reconstructedStream = [];
% conditionCounts = zeros(size(reconstructTrials,2), size(reconstructTrials,3));
% for i = 1:size(rawTrials,3)
%     conditionID = mod((i-1),6)+1;
%     
%     F1ID = mod(conditionID - 1,3)+1;
%     F2ID = floor((conditionID-1)/3)+1;
%     
% %     thisIDs = (i-1)*trialLength + (1:trialLength);
% %     thisLoopIDs = app.SavedData.BestStateMap(thisIDs,1);
% %     thisPhaseIDs = app.SavedData.BestStateMap(thisIDs,2);
%     
% %     thisStateIDs = [];
% %     for j = 1:length(thisLoopIDs)
% %         thisStateIDs(j) = stateLookup(thisLoopIDs(j), thisPhaseIDs(j));
% %     end
%     
%     conditionCounts(F1ID, F2ID) = conditionCounts(F1ID, F2ID) + 1;
%     rawConditions(F1ID, F2ID,:,:,conditionCounts(F1ID, F2ID)) = rawTrials(:,:,i);
%         
%     for j = 1:10
%         currentCluster = allClusters(i);
% 
%         for t = 2:simulateTime
%             currentCluster(t) = randsample(1:size(bestModel,1), 1, true, bestModel(currentCluster(t-1),:));
%         end
% 
%         reconstructedStream(:, F1ID, F2ID, :, conditionCounts(F1ID, F2ID), j) = clusterMeans(currentCluster,:)';
%     end
% end
% 
% %% Calc correlation
% 
% % Rsquareds = [];
% % for i = 1:size(reconstructedStream,2)
% %     for j = 1:size(reconstructedStream,3)
% %         for k = 1:size(reconstructedStream,1)
% %             for trialNum = 1:size(reconstructedStream,5)
% %                 thisReconstruction = squeeze(reconstructedStream(k,i,j,:,trialNum));
% %                 targetData = squeeze(reconstructTrials(k,i,j,:));
% %                 
% %                 Rsquareds(k,i,j,trialNum) = corr(thisReconstruction, targetData);
% %             end
% %         end
% %     end
% % end
% % 
% % meanRsquared = nanmean(Rsquareds(:))
% 
% useTimes = startIndex+1:size(rawConditions,4);
% 
% trialRsquareds = [];
% for i = 1:size(reconstructedStream,2)
%     for j = 1:size(reconstructedStream,3)
%         for trialNum = 1:size(reconstructedStream,5)
%             targetData = squeeze(rawConditions(i,j,:,useTimes,trialNum));
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
%             targetData = squeeze(reconstructTrials(:,i,j,useTimes));
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
% 
% F1Start = 0;
% F2Start = 3.5;
% trialEnd = 5.5;
% plotF1 = 1;
% plotF2 = 1;
% plotTrial = 3;
% trialEnd = 5.5;
% 
% startTimes = 1:startIndex;
% 
% thisStart = squeeze(reconstructTrials(:,plotF1,plotF2,startTimes));
% thisTarget = squeeze(reconstructTrials(:,plotF1,plotF2,useTimes));
% thisSimulated = squeeze(reconstructedStream(:,plotF1,plotF2,:,plotTrial, :));
% 
% totalFiringRates = abs(pcaBasis(:,1));
% [~, topFiring] = sort(totalFiringRates, 'descend');
% 
% plotNeurons = topFiring([6,15,21]);
% 
% originalTime = time;
% originalTime = decimate(originalTime(times),decimateAmount);
% originalTime = originalTime(end-1-trialLength:end-2);
% 
% colors = lines(100);
% 
% figureHandle = figure(8);
% figureHandle.Renderer='Painters';
% clf;
% hold on;
% for i = 1:length(plotNeurons)
%     plot(originalTime([startTimes useTimes(1)]), [thisStart(plotNeurons(i),:) thisTarget(plotNeurons(i),1)], 'Color', 'k');
% end
% 
% for i = 1:length(plotNeurons)
%     plot(originalTime(useTimes), thisTarget(plotNeurons(i),:), 'Color', colors(i,:), 'LineWidth', 3);
% end
% 
% for i = 1:length(plotNeurons)
%     for j = 1:size(thisSimulated,3)
%         plot(originalTime(useTimes), thisSimulated(plotNeurons(i),:,j), 'Color', colors(i,:));
%     end
% end
% xlim([min(originalTime), trialEnd])


