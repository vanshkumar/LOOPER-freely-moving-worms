function [correlation, constraints] = ComapreToMarkovFunctionLOOPER(nn, repop, stateCounts)

   try       
    %% Load validation data FOR MONKEY
    
    currentTime = clock;
    rng(floor(currentTime(6)*1000));

    clear time;

    romoAnalysisPipeline_forFig3
    forValidation = 1;
    monkeySetupForFig3
    % load('monkeyDataForFigure3')
    
    load('monkeyDataForFigure3.mat');

    app.SavedData = saveData;
    
    originalData = saveData;

    numTrial = max(saveData.TrialData);
    originalLength = size(saveData.RawData,2)/numTrial;
    trialLength = size(saveData.FinalStream,1) / numTrial;
    
    rawTrials = [];
    for i = 1:size(finalTrials,3)
        rawTrials(:,:,i) = dataPCABAsis * finalTrials(:,:,i);
    end

    trialData = convertToCell(rawTrials(:,1+originalLength-trialLength:originalLength,:)) ;  
    time = timeInt(1)/1000 : 0.01 : timeInt(2)/1000;
    times = 1:900;
    
    originalTime = time(times);
    originalTime = decimate(originalTime,decimateAmount);
    originalTime = originalTime(end-1-trialLength:end-2);

    startTime = 4.0;
    [~, startIndex] = min(abs(originalTime - startTime));

    F1Start = 0;
    F2Start = 3.5;
    trialEnd = 5.5;
    plotTrial = 50;
    trialEnd = 5.5;
    
    %% Load LOOPER data

    numTrial = max(app.SavedData.TrialData);
    trialLength = size(app.SavedData.FinalStream,1) / numTrial;
    trialStarts = 1:trialLength:size(app.SavedData.FinalStream,1);
    
    params = [];
    params.NearestNeighbors = nn;
    params.RepopulateDensity = repop;
    params.TotalStates = stateCounts;
    
    LOOPER(originalData, true, [], [], [], params);  

    allLOOPERStates = [];
    for i = 1:size(app.SavedData.BestStateMap,1)
        allLOOPERStates(i) = find(ismember(app.SavedData.BestLoopAssignments, app.SavedData.BestStateMap(i,:), 'rows'));
    end

    %% Simulate traces

    simulateTime = size(finalTrials,2) - startIndex;

    bestStates = allLOOPERStates;
    bestModel = app.SavedData.BestModel;

    finalIDs = [];
    for i = 1:length(trialData)
        thisIDs = (i-1) * originalLength + (1+originalLength-trialLength:originalLength);
        finalIDs = [finalIDs thisIDs];
    end

    finalStream = app.SavedData.RawData(:,finalIDs)' * dataPCABAsis';

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
        for j = 1:size(app.SavedData.BestEmission,1)
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

    trialRsquareds = [];
    for i = 1:size(reconstructedStream,2)
        targetData = squeeze(trialData{i}(:,useTimes));

        for repeatCount = 1:size(reconstructedStream,4)
            thisReconstruction = squeeze(reconstructedStream(:,i,1:SIM_TIME, repeatCount));


            trialRsquareds(i,repeatCount) = corr(thisReconstruction(:), targetData(:));
        end
    end


    meanTrialRsquared = nanmean(trialRsquareds(:))
    
    
    correlation = 1 - abs(meanTrialRsquared);
    constraints = -1;
   catch ME
    warning(ME.message);
    disp([ME.stack(1).name ' Line: ' num2str(ME.stack(1).line)]);
    correlation = 1;
    constraints = 1;
   end
end

