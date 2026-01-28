currentTime = clock;
rng(floor(currentTime(6)*1000));

allTrials = {};
trialIndex = 1;
singleTrace = [];
for i = 1:2
    for j = 1:totalTrials

        if ~SHOULD_VALIDATE
            thisIDs = useIDs;
        else
            thisIDs = validationIDs;
        end

        if i == 1
            movieTarget = 2;
            movieOfftarget = 4;
        else
            movieTarget = 4;
            movieOfftarget = 2;
        end

        moviePoolCorrect = find(ismember(eventData.subjID, thisIDs) & ...
            (cell2mat(eventData.resp) == movieTarget & cell2mat(eventData.respAcc) == 0) & ...
            eventData.movement == 0 & contains(eventData.matfile, 'wholeBrainPower'));
        moviePoolWrong = find(ismember(eventData.subjID, thisIDs) & ...
            (cell2mat(eventData.resp) == movieOfftarget & cell2mat(eventData.respAcc) == 2) & ...
            eventData.movement == 0 & contains(eventData.matfile, 'wholeBrainPower'));

        numIncorrectBootstraps = floor(numBootstraps * INCORRECT_RATIO);
        numCorrectBootstraps = numBootstraps - numIncorrectBootstraps;

        correctBootstraps = randsample(length(moviePoolCorrect), numCorrectBootstraps, true);
        incorrectBootstraps = randsample(length(moviePoolWrong), numIncorrectBootstraps, true);

        bootstrapIDs = [moviePoolCorrect(correctBootstraps); moviePoolWrong(incorrectBootstraps)];

        thisData = [];
        for k = 1:length(bootstrapIDs)
%             files = dir([loadPath eventData.matfile{bootstrapIDs(k)}]);
% 
%             load([loadPath files(1).name]);

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