%%
% eventData = csv2struct('../../6_1200subj_eventFiles\subj_onset_behav_1200subjs_0413.csv');
eventData = csv2struct('\\192.168.1.206\LabData\Connor\HCPData\Social_cog_task\6_1200subj_eventFiles\subj_onset_behav_1200subjs_0413_withMovemet.csv');
loadPath = '\\192.168.1.206\LabData\Connor\HCPData\Social_cog_task\5_1200subjs_time_series\power_2011_atlas_mat\';

%%

waitHandle = parfor_progressbar(length(eventData.matfile), 'Loading all data');


preloadedData = {};
loadedID = 1;
for i = 1:length(eventData.matfile)
    
    files = dir([loadPath eventData.matfile{i}]);

    waitHandle.iterate(1);
    
    if isempty(files)
        continue;
    end
    
    load([loadPath files(1).name]);
    
    preloadedData{loadedID} = arr;
    
    eventData.preloadID(i) = loadedID;
    loadedID = loadedID + 1;
end
waitHandle.close();

%%

SHOULD_VALIDATE = false;
VALDIATION_AMOUNT = 0.3;
INCORRECT_RATIO = 0.3;
numBootstraps = 150;
totalTrials = 10;
PRE_MOVIE_FRAMES = 0;
POST_MOVIE_FRAMES = 0;

movies = unique(eventData.movie);
subjects = unique(eventData.subjID);

subjects(subjects == 171734) = [];
subjects(subjects == 103010) = [];
subjects(subjects == 748662) = [];
subjects(subjects == 175540) = [];



clusterCounts = [16, 20, 24];
delayLengths = [2, 3, 4];
smoothings = [1.5, 2, 2.5];
nearestNeighborCounts = [4, 5, 6];
totalTrialCounts = 10;
totalValidationCounts = 10;

dataLength = 25 + PRE_MOVIE_FRAMES + POST_MOVIE_FRAMES;

totalParameterTests = length(clusterCounts) * length(delayLengths) * length(smoothings) * length(nearestNeighborCounts) * totalTrialCounts;

allLoopCounts = nan(length(clusterCounts), length(delayLengths), length(smoothings), length(nearestNeighborCounts), totalTrialCounts, dataLength);
allLoopAccuracies = nan(length(clusterCounts), length(delayLengths), length(smoothings), length(nearestNeighborCounts), totalTrialCounts, totalValidationCounts, dataLength);
allMeanAccuracies = nan(length(clusterCounts), length(delayLengths), length(smoothings), length(nearestNeighborCounts), totalTrialCounts, dataLength);

counter = 0;
for clusterCountID = 1:length(clusterCounts)
    for delayLengthID = 1:length(delayLengths)
        for smoothingID = 1:length(smoothings)
            for nearestNeighborsID = 1:length(nearestNeighborCounts)
                for dataTrial = 1:totalTrialCounts
                    
                    rng('default');
%                     rng(1);
                    currentTime = clock;
                    rng(floor(currentTime(6)*1000));

                    allIDs = subjects;
                    useIDs = randsample(subjects, ceil(length(subjects)*VALDIATION_AMOUNT), false);
                    validationIDs = setdiff(allIDs, useIDs);

                    
                    clusterCount = clusterCounts(clusterCountID);
                    delayLength = delayLengths(delayLengthID);
                    smoothing = smoothings(smoothingID);
                    nearestNeighbors = nearestNeighborCounts(nearestNeighborsID);
                    
                    disp(['Building data...']);
                    SHOULD_VALIDATE = false;
                    
                    getData;


                    load('F:\Dropbox\MethodsPaperCode\Data\socialMovies2.mat');

                    params = [];
                	params.PreprocessData.DelayCount = delayLength;
                	params.PreprocessData.Smoothing = smoothing;
                    params.NearestNeighbors = nearestNeighbors;
                    params.PutativeClusterCounts = [clusterCount];
                    params.PutativeLoopCounts = [2];
                    params.RepopulateDensity = 0.5;

                    LOOPER(saveData, true, allTrials, [], [], params);
                    
                    validateData;
                    
                    thisLoopCounts = max(loopAssignments);
                    if range(thisLoopCounts) == 0
                        thisLoopCounts = min(loopAssignments);
                        thisLoopCounts = 3 - thisLoopCounts;
                    end
                    
                    thisDataLength = min(length(thisLoopCounts), dataLength);
                    
                    allLoopCounts(clusterCountID, delayLengthID, smoothingID, nearestNeighborsID, dataTrial, end-thisDataLength+1:end) = thisLoopCounts(end-thisDataLength+1:end);
                    modelLoopAssignments = loopAssignments;
                    
                    %%
                    
                    disp('Validating data...');
                    
                    SHOULD_VALIDATE = true;
                    
                    thisModelAccuracies = nan(totalValidationCounts, thisDataLength);
                    for valdiationID = 1:totalValidationCounts
                        getData;

                        validateData;
                        
                        allAccuracies = modelLoopAssignments == loopAssignments;
                        thisAccuracies = mean(allAccuracies);
                        
                        thisModelAccuracies(valdiationID, end-thisDataLength+1:end) = thisAccuracies(end-thisDataLength+1:end);
                        allLoopAccuracies(clusterCountID, delayLengthID, smoothingID, nearestNeighborsID, dataTrial, valdiationID, end-thisDataLength+1:end) = thisAccuracies(end-thisDataLength+1:end);
                    end
                    
                    allMeanAccuracies(clusterCountID, delayLengthID, smoothingID, nearestNeighborsID, dataTrial, end-thisDataLength+1:end) = mean(thisModelAccuracies);
                    
                    counter = counter + 1;
                    disp(['Finsihed ' num2str(counter) ' of ' num2str(totalParameterTests)]);
                    
                end
            end
        end
    end
end

%%






allWeights = allMeanAccuracies/nansum(allMeanAccuracies(:))
weightedLoopCounts = allLoopCounts .* allWeights;

weightedLoopCounts = nansum(reshape(weightedLoopCounts, [], size(weightedLoopCounts,6))) * 25;

plot(weightedLoopCounts)

