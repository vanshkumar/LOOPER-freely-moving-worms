
%% Load validation data FOR MONKEY

load('monkeyDataForFigure3.mat');

app.SavedData = saveData;


romoAnalysisPipeline_forFig3
forValidation = 1;
NUM_BOOTSTRAPS = 10;
monkeySetupForFig3
% load('monkeyDataForFigure3')

numTrial = max(saveData.TrialData);
originalLength = size(saveData.RawData,2)/numTrial;
trialLength = size(saveData.FinalStream,1) / numTrial;

rawTrials = [];
for i = 1:size(finalTrials,3)
    rawTrials(:,:,i) = dataPCABAsis * finalTrials(:,:,i);
end

trialData = convertToCell(rawTrials(:,1+originalLength-trialLength:originalLength,:)) ;

rawData = [];
for i = 1:size(finalTrials,3)
    rawData(:,:,i) = finalTrials(:,:,i);
end

rawData = convertToCell(rawData(:,3:originalLength,:)) ;

lastEnd = 0;
rawTrialData = [];
for i = 1:length(rawData)
    rawTrialData(lastEnd + (1:size(rawData{i}, 2))) = i;

    lastEnd = lastEnd + size(rawData{i}, 2);
end

rawData = mergeData(rawData);

[processedData, processedTrialData, procesedTrialSwitches] = preprocessData(rawData, [], 0, [], 0, rawTrialData, 0, true, saveData.PreprocessData.Smoothing, saveData.PreprocessData.ZScore, saveData.PreprocessData.DelayTime, saveData.PreprocessData.DelayCount, saveData.DataMean, saveData.DataSTD);


originalTime = time;
originalTime = decimate(originalTime(times),decimateAmount);
originalTime = originalTime(end-1-trialLength:end-2);

startTime = 4.0;
[~, startIndex] = min(abs(originalTime - startTime));

F1Start = 0;
F2Start = 3.5;
trialEnd = 5.5;
plotTrial = 50;
trialEnd = 5.5;

%%

load('monkeyDataForFigure3.mat');

app.SavedData = saveData;

SHOULD_VALIDATE = 0;
Figure3PlotLoops

load('BayesoptResultsLOOPER.mat')

nn = BayesoptResults.XAtMinObjective.nn;
repop = BayesoptResults.XAtMinObjective.repop;
stateCounts = BayesoptResults.XAtMinObjective.stateCounts;

numTrial = max(app.SavedData.TrialData);
trialLength = size(app.SavedData.FinalStream,1) / numTrial;
trialStarts = 1:trialLength:size(app.SavedData.FinalStream,1);

originalData = app.SavedData;

trialConditions = mod(0:numTrial-1, 6)+1;


%% Build LOOPER

params = [];
params.NearestNeighbors = nn;
params.RepopulateDensity = repop;
params.TotalStates = stateCounts;

LOOPER(originalData, true, [], [], [], params); 

app.SavedData = saveData;

%%

allLOOPERStates = [];
for i = 1:size(originalData.BestStateMap,1)
    allLOOPERStates(i) = find(ismember(originalData.BestLoopAssignments, originalData.BestStateMap(i,:), 'rows'));
end

allHMMStates = [];
for i = 1:size(app.SavedData.BestStateMap,1)
    allHMMStates(i) = find(ismember(app.SavedData.BestLoopAssignments, app.SavedData.BestStateMap(i,:), 'rows'));
end


%% Reconstruction

for USE_HMM = 0:1
    if USE_HMM
        bestStates = allHMMStates;
        bestModel = app.SavedData.BestModel;
    else
        bestStates = allLOOPERStates;
        bestModel = originalData.BestModel;
    end

    finalIDs = [];
    for i = 1:length(trialData)
        thisIDs = (i-1) * originalLength + (1+originalLength-trialLength:originalLength);
        finalIDs = [finalIDs thisIDs];
    end

    finalStream = originalData.RawData(:,finalIDs)' * dataPCABAsis';

    clusterMeans = [];
    for i = 1:size(bestModel,1)
        thisIndicies = find(bestStates == i);

        if isempty(thisIndicies)
            clusterMeans(i,:) = zeros(1,size(finalStream,2));
        else
            clusterMeans(i,:) = mean(finalStream(thisIndicies,:), 1);
        end
    end

    reconstructedStream = [];
    for i = 1:length(trialData)
        thisTrialIDs = (i-1)*trialLength + (1:trialLength);
        reconstructedStream(:, i, :) = clusterMeans(bestStates(thisTrialIDs),:)';
    end

    %% Calc correlation

    trialRsquareds = [];
    for i = 1:size(reconstructedStream,2)
        thisTrialIDs = (i-1)*trialLength + (1:trialLength);
        targetData = finalStream(thisTrialIDs,:)';

        thisReconstruction = squeeze(reconstructedStream(:,i, :));

        trialRsquareds(i) = corr(thisReconstruction(:), targetData(:));
    end


    meanTrialRsquared = nanmedian(trialRsquareds(:))

    figureHandle = figure(9);
    figureHandle.Renderer='Painters';
    clf;
    ksdensity(trialRsquareds(:));

    if USE_HMM
        hmmDistribution = trialRsquareds;
    else
        looperDistribution =  trialRsquareds;
    end
end

%%

if exist('hmmDistribution') && exist('looperDistribution')
    figure(1);
    clf;
    hold on;
    ksdensity(hmmDistribution(:));
    ksdensity(looperDistribution(:));
    
    p = ranksum(mean(hmmDistribution,2), mean(looperDistribution,2))
    
    looperR2 = median(looperDistribution);
%     [~,looperNumPCs] = min(abs(cumsum(dataExplained) - looperR2*100));
end


%% Simulate traces

for USE_HMM = 0:1
    simulateTime = size(finalTrials,2) - startIndex;

    if USE_HMM
        bestStates = allHMMStates;
        bestModel = app.SavedData.BestModel;
        % clusterMeans = ((means' .* streamSTD) + streamMean)' * pcaBasis';
    else
        bestStates = allLOOPERStates;
        bestModel = originalData.BestModel;
        % clusterMeans = app.SavedData.BestEmission;
    end

    finalIDs = [];
    for i = 1:length(trialData)
        thisIDs = (i-1) * originalLength + (1+originalLength-trialLength:originalLength);
        finalIDs = [finalIDs thisIDs];
    end

    finalStream = originalData.RawData(:,finalIDs)' * dataPCABAsis';

    clusterMeans = [];
    for i = 1:size(bestModel,1)
        thisIndicies = find(bestStates == i);

        if isempty(thisIndicies)
            clusterMeans(i,:) = zeros(1,size(finalStream,2));
        else
            clusterMeans(i,:) = mean(finalStream(thisIndicies,:), 1);
        end
    end

    allClusters = [];
    for i = 1:length(trialData)
        thisStart = trialData{i}(:,startIndex);

        clusterDistances = [];
        for j = 1:size(originalData.BestEmission,1)
            thisIndicies = find(bestStates == j);

            stds = std(finalStream(thisIndicies,:), [], 1);        
            thisEmission = mean(finalStream(thisIndicies,:), 1);

            clusterStream = (thisEmission - thisStart') ./ stds;           

            clusterDistances(j,:) = sum(clusterStream.^2,2);
        end

        [~, bestClusters] = min(clusterDistances, [], 1);

        allClusters(i) = bestClusters;        

    end

    reconstructedStream = [];
    for i = 1:length(allClusters)            
        for j = 1:50
            currentCluster = allClusters(i);

            for t = 2:simulateTime
                currentCluster(t) = randsample(1:size(bestModel,1), 1, true, bestModel(currentCluster(t-1),:));
            end

            reconstructedStream(:, i, :, j) = clusterMeans(currentCluster,:)';
        end
    end

    %% Calc correlation

    SIM_TIME = 19;
    useTimes = startIndex:startIndex+SIM_TIME-1;

    % figure(2);
    % clf;
    % hold on;
    % for i = 1:size(reconstructedStream,2)
    %     thisTrace = squeeze(reconstructedStream(:,i,:,1));
    %     
    %     targetData = squeeze(trialData{i});
    %     
    % %     plot3(thisTrace(1,:), thisTrace(2,:), thisTrace(3,:));
    %     plot3(targetData(1,:), targetData(2,:), targetData(3,:));
    % end

    trialRsquareds = [];
    for i = 1:size(reconstructedStream,2)
        targetData = squeeze(trialData{i}(:,useTimes));

        for repeatCount = 1:size(reconstructedStream,4)
            thisReconstruction = squeeze(reconstructedStream(:,i,1:SIM_TIME, repeatCount));

            trialRsquareds(i,repeatCount) = corr(thisReconstruction(:), targetData(:));
        end
    end


    meanTrialRsquared = nanmedian(trialRsquareds(:))

    figureHandle = figure(9);
    figureHandle.Renderer='Painters';
    clf;
    ksdensity(trialRsquareds(:));

    if USE_HMM
        hmmDistribution = trialRsquareds;
    else
        looperDistribution =  trialRsquareds;
    end
end

%%

if exist('hmmDistribution') && exist('looperDistribution')
    figure(1);
    clf;
    hold on;
    ksdensity(hmmDistribution(:));
    ksdensity(looperDistribution(:));
    
    p = ranksum(mean(hmmDistribution,2), mean(looperDistribution,2))
    
%     ps = [];
%     for i = 1:size(hmmDistribution,1)
%         ps(i) = ranksum(hmmDistribution(i,:), looperDistribution(i,:));
%     end
end

%% Decode data

decodeRateLOOPER = [];
decodeRateHMM = [];
decodeRatio = [];

for bootstrapTrial = 1:10
    for USE_HMM = 0:1
        if USE_HMM
            bestStates = allHMMStates;
        else
            bestStates = allLOOPERStates;
        end

        conditionOccupancies = zeros(max(bestStates), trialLength, 6);
        allStates = [];
        for i = 1:length(trialData)
            sameConditions = find(trialConditions == trialConditions(i));
            bootstrapID = sameConditions(floor(rand(1)*length(sameConditions)) + 1);

            thisIndicies = (bootstrapID-1)*trialLength + (1:trialLength);

            bestSequence = bestStates(thisIndicies);

            allStates(i,:) = bestSequence;

            for j = 1:length(bestSequence)
                thisState = bestSequence(j);
                conditionOccupancies(thisState, j, trialConditions(bootstrapID)) = conditionOccupancies(thisState, j, trialConditions(bootstrapID)) + 1;
            end
        end

        conditionOccupanciesStateMarginalized = squeeze(sum(conditionOccupancies,2));
        conditionOccupanciesStateMarginalized = conditionOccupanciesStateMarginalized ./ sum(conditionOccupanciesStateMarginalized,2);
        conditionOccupanciesStateMarginalized(isnan(conditionOccupanciesStateMarginalized)) = 0;


        %% Esitamte decoding rate
        times = 1:900;

        originalTime = time;
        originalTime = decimate(originalTime(times),decimateAmount);
        originalTime = originalTime(end-1-trialLength:end-2);

        truePlotTime = originalTime;

        bestLoopCounts = [];
        for i = 1:size(trialLoopIDs,1)    
            bestLoopCounts(i) = length(unique(trialLoopIDs(i,:)));
        end

        F1StartTime = truePlotTime(find(bestLoopCounts == 3, 1));
        F1EndTime = 3;
        F2StartTime = truePlotTime(find(bestLoopCounts == 6, 1));
        F2EndTime = truePlotTime(find(bestLoopCounts(1:end-1) == 6 & bestLoopCounts(2:end) ~= 6, 1));



        loopOccupancies = conditionOccupancies ./ sum(conditionOccupancies,3);
        loopOccupancies(isnan(loopOccupancies)) = 0;

        F1Indices = find(truePlotTime >= F1StartTime & truePlotTime <= F1EndTime);
        F2Indices = find(truePlotTime >= F2StartTime & truePlotTime <= F2EndTime);

        F1Indices(F1Indices > trialLength) = [];
        F2Indices(F2Indices > trialLength) = [];

        F1Rates = [];
        F2Rates = [];
        for i = 1:6
            conditionID = i;

            totalDecodeRates = [];
            for j = 1:trialLength
                loopDistribution = conditionOccupancies(:,j,i);
                loopDistribution = loopDistribution / sum(loopDistribution);

                totalDecodeRates(j) = sum(loopDistribution .* loopOccupancies(:,j,i));                
            end

            F1Rates(i) = mean(abs(totalDecodeRates(F1Indices) - 0.5))*2;
            F2Rates(i) = mean(abs(totalDecodeRates(F2Indices) - 1));
        end
        F1Rate = 1 - mean(F1Rates);
        F2Rate = 1 - mean(F2Rates);

        if USE_HMM == 1
            decodeRateHMM(bootstrapTrial) = (F1Rate + F2Rate) / 2;
        else
            decodeRateLOOPER(bootstrapTrial) = (F1Rate + F2Rate) / 2;
        end
            

        %% Calculate P(condition | time, true condition)

        allTimeProbabilities = [];
        for i = 1:trialLength
            timeStates = allStates(:,i);

            for j = 1:length(unique(trialConditions))
                thisConditions = find(trialConditions == j);

                thisTimeProbabilities = conditionOccupanciesStateMarginalized(timeStates(thisConditions),:);
                thisTimeProbabilities = nanmean(thisTimeProbabilities);

                allTimeProbabilities(j, i,:) = thisTimeProbabilities;
            end
        end

        % This is P(CORECT condition | time, true condition)
        figure(10 + USE_HMM);
        clf;
        hold on;
        for j = 1:length(unique(trialConditions))
            plot(truePlotTime, squeeze(allTimeProbabilities(j,:,j)));
        end

        plot(ones(2,1)*F1Start, [0 1], 'r');
        plot(ones(2,1)*F2Start, [0 1], 'r');
    end
    
    decodeRatio(bootstrapTrial) = (decodeRateLOOPER - 0.25) / (decodeRateHMM - 0.25);
end

%%

labels = {'LOOPER Simulation R^2', 'LOOPER BayesOpt Simulation R^2', 'LOOPER Decode Rate', 'LOOPER BayesOpt Decode Rate'};
ticks = [1 2 3 4];

figure(1);
clf;
hold on;

% meanValues = [mean(looperDistribution(:)) mean(hmmDistribution(:)) mean(decodeRateLOOPER(:)) mean(decodeRateHMM)];
% errorLow = [-std(looperDistribution(:))/sqrt(length(looperDistribution(:)))*1.98 -std(hmmDistribution(:))/sqrt(length(hmmDistribution(:)))*1.98 -std(decodeRateLOOPER), -std(decodeRateHMM) ];
% errorHigh = [std(looperDistribution(:))/sqrt(length(looperDistribution(:)))*1.98 std(hmmDistribution(:))/sqrt(length(hmmDistribution(:)))*1.98 std(decodeRateLOOPER),  std(decodeRateHMM) ];
% 
% bar(ticks,meanValues)
% er = errorbar(ticks,meanValues,errorLow,errorHigh); 
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  

boxplot([mean(looperDistribution,2)' mean(hmmDistribution,2)' decodeRateLOOPER decodeRateHMM], [ones(1, length(mean(looperDistribution,2)))*ticks(1) ones(1, length( mean(hmmDistribution,2)))*ticks(2) ones(1, length(decodeRateLOOPER(:)))*ticks(3) ones(1, length(decodeRateHMM(:)))*ticks(4)]);

set(gca,'xtick',ticks)
set(gca,'xticklabel',labels)
xtickangle(45)
