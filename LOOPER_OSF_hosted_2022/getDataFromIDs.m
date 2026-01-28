currentTime = clock;
rng(floor(currentTime(6)*1000));

allTrials = {};
trialIndex = 1;
singleTrace = [];
for i = 1:length(conditionIDs)
    for j = 1:totalTrials

        thisConditionIDs = conditionIDs{i};

        bootstraps = randsample(length(thisConditionIDs), numBootstraps, true);

        bootstrapIDs = thisConditionIDs(bootstraps);

        thisData = [];
        for k = 1:length(bootstrapIDs)
            arr = preloadedData{eventData.preloadID(bootstrapIDs(k))};

            startTime = round(eventData.taskOnsetTR(bootstrapIDs(k))) - PRE_MOVIE_FRAMES;
            endTime = startTime + round(20/0.72) + PRE_MOVIE_FRAMES + POST_MOVIE_FRAMES;

            thisData(:,:,k) = arr(startTime:endTime,:);
        end
        thisData = mean(thisData,3);

        allTrials{trialIndex} = (thisData)';
        singleTrace = [singleTrace thisData'];

        trialIndex = trialIndex + 1;
    end
end