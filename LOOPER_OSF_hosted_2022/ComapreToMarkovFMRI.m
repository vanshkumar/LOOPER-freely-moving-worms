%% Load validation data FOR FMRI

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
    
%% Run hmm fits
    
% This data was bult using USE_VALIDATION = true, so it needs to be false
% for validation below
load('socialMovieSingle.mat');

originalSaveData = saveData;

app.SavedData = originalSaveData;

numTrial = max(app.SavedData.TrialData);
trialLength = size(app.SavedData.FinalStream,1) / numTrial;
trialStarts = 1:trialLength:size(app.SavedData.FinalStream,1);
% trialStarts = [trialStarts size(app.SavedData.FinalStream,1)+1];

%pcs=11 eps=4 t=11 k=60 #=105 Bad?
%pcs=9 eps=5 t=13 k=60 #=105 Less bad?
%pcs=14 eps=5.3477 t=1 k=18 #=100 pattern search

load('BayesoptResultsFMRI.mat')

pcCount = BayesoptResults.XAtMinObjective.pcCount;
epsilon = BayesoptResults.XAtMinObjective.epsilon;
tValue = BayesoptResults.XAtMinObjective.t;
kValue = BayesoptResults.XAtMinObjective.k;
numStates = BayesoptResults.XAtMinObjective.stateCounts;

allFits = [];
bestFit = 0;
bestFitID = 0;
allHMMHistograms = [];
allLOOPERHistograms = [];
epsilons = ones(10,1) * epsilon;
ts = ones(10,1) * tValue;
ks = ones(10,1) * kValue;
Qvalues = ones(10,1) * numStates;
pcCounts = ones(10,1) * pcCount;

for parameterSets = 1%:length(Qvalues)
    
    dimensionalityStream = app.SavedData.FinalStream;

    pcaBasis = eye(size(dimensionalityStream,2));
    
    [pcaBasis, ~, ~, ~, explained] = pca(dimensionalityStream);
%     pcaCutoff = find(cumsum(explained) > explainedVariances(parameterSets), 1);
    pcaCutoff = pcCounts(parameterSets);
    pcaBasis = pcaBasis(:, 1:pcaCutoff);

    dimensionalityStream = dimensionalityStream * pcaBasis;
    
    rawStream = dimensionalityStream;

    [rawStream, diffusionMap] = diffusionMapReduction(dimensionalityStream, epsilons(parameterSets), ts(parameterSets), ks(parameterSets));

%     startData = tempData' * pcaBasis;
%     
%     [~, newStream] = diffusionMapTransformation(startData, dimensionalityStream, diffusionMap, epsilons(parameterSets), ts(parameterSets), ks(parameterSets));
    
%     figure(1);
%     clf;
%     hold on;
%     plot3(startData(:,1),startData(:,2),startData(:,3));
%     plot3(newStream(:,1),newStream(:,2),newStream(:,3));
    
%     [rawStream, diffusionMap] = diffusionMapTransformation(dimensionalityStream, epsilons(parameterSets), ts(parameterSets), ks(parameterSets));
%             [~, originalData] = diffusionMapTransformation(((means' .* streamSTD) + streamMean)', dimensionalityStream, diffusionMap, epsilons(parameterSets), ts(parameterSets), ks(parameterSets));
%             clusterMeans = originalData * pcaBasis';
%     rawStream = dimensionalityStream;
%     rawStream = diffusion_maps(dimensionalityStream, 10, 1/2, 0.1);
    
    

    streamMean = mean(rawStream', 2);
    streamSTD = std(rawStream', [], 2);

    thisIndex = 1;
    data = [];
    for i = 1:numTrial
        finalindices = thisIndex:thisIndex+trialLength-1;
        thisIndex = thisIndex + trialLength;

        data(:,:,i) = (rawStream(finalindices,1:size(rawStream,2))' - streamMean) ./ streamSTD;
    end
    
    O = size(data,1);          %Number of coefficients in a vector 
    T = size(data,2);         %Number of vectors in a sequence 
    nex = size(data,3);        %Number of sequences 
    M = 1;          %Number of mixtures 
    Q = Qvalues(parameterSets);%floor(size(app.SavedData.BestModel,1) / 1.2^(parameterSets-1));          %Number of states 
    cov_type = 'diag';
    
    %% Build LOOPER
    
    % nn = 
    
%     nn = 11;
%     repopulationDensity = 0.99;
%     stateCounts = 203;
%     
%     params = [];
%     params.NearestNeighbors = nn;
%     params.RepopulateDensity = repopulationDensity;
%     params.TotalStates = stateCounts;
%     
%     LOOPER(originalSaveData, true, [], [], [], params);
%     LOOPER(originalSaveData, true);
%     
%     app.SavedData = saveData;

    %% Build HMM

    allPercents = [];
    CPT1s = [];
    CPT3s = [];
    allMeans = [];
    for runNumber = 1:1

        includeRows = 1;
        
        matrixPrior = zeros(Q,Q);
        if includeRows > 0
            matrixPrior = diag(ones(Q-1,1),1) * includeRows;
            matrixPrior(Q,1) = includeRows;
            
            for i = 2:includeRows
                matrixPrior = matrixPrior + diag(ones(Q-i,1),i) * (includeRows - i);
%             matrixPrior(Q,i) = (includeRows - i);
            end
        end
%         matrixPrior = kron(eye(3),matrixPrior) * 5000/Q/3;
%         matrixPrior = matrixPrior + eye(size(matrixPrior)) * 1;


        intra = zeros(2);
        intra(1,2) = 1; % node 1 in slice t connects to node 2 in slice t

        inter = zeros(2);
        inter(1,1) = 1; % node 1 in slice t-1 connects to node 1 in slice t

        ns = [Q O];
        dnodes = [1];
        onodes = [2];

        eclass1 = [1 2];
        eclass2 = [3 2];
        eclass = [eclass1 eclass2];

        covariates = zeros(O,O,Q);
        for i = 1:Q
            covariates(:,:,i) = eye(O)*0.2;
        end

        bnet = mk_dbn(intra, inter, ns, 'discrete', dnodes, 'observed', onodes, 'eclass1', eclass1, 'eclass2', eclass2);
        prior0 = normalise(ones(Q,1));
        % transmat0 = app.SavedData.BestModel;
        transmat0 = mk_stochastic(rand(Q,Q));
        obsmat0 = mk_stochastic(rand(Q,O));
        % obsmat0 = ((app.SavedData.BestEmission' - streamMean) ./ streamSTD)';
        obsmat0(isnan(obsmat0)) = 0;
        bnet.CPD{1} = tabular_CPD(bnet, 1, prior0);
        bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', obsmat0, 'cov', covariates, 'cov_type', 'diag', 'clamp_cov', true);
        bnet.CPD{3} = tabular_CPD(bnet, 3, 'CPT', transmat0, 'prior_type', 'dirichlet', 'dirichlet_type', ...
           'prior', 'dirichlet_weight', matrixPrior);
        % bnet.CPD{3} = tabular_CPD(bnet, 3, 'CPT', transmat0, 'prior_type', 'dirichlet', 'dirichlet_type', ...
        %    'prior', 'dirichlet_weight', matrixPrior);

        engine = bk_inf_engine(bnet, 'clusters', 'exact');
        % engine = hmm_inf_engine(bnet);

        ss = 2;%slice size(ss)
        max_iter=500;%iterations for EM
        cases = cell(1, nex);
        for i=1:nex
          cases{i} = cell(ss,T);
        %   cases{i}(onodes,:) = ev(onodes, :);
          for t = 1:T
            cases{i}{2,t} = squeeze(data(:, t, i));
          end
        end

        [bnet2, LLtrace] = learn_params_dbn_em(engine, cases, 'max_iter', max_iter);

        CPT1 = get_field(bnet2.CPD{1}, 'cpt');
        means = get_field(bnet2.CPD{2}, 'mean')';
        % means = (pcaBasis * ((get_field(bnet2.CPD{2}, 'mean') .* streamSTD) + streamMean))';
        covs = get_field(bnet2.CPD{2}, 'cov');
        CPT3 = get_field(bnet2.CPD{3}, 'cpt');

%         figure(1);
%         clf;
%         subplot(2,1,1);
%         imagesc(means);
%         subplot(2,1,2);
%         imagesc(app.SavedData.BestEmission);


        %% Compare models
        % CPT1 = get_field(bnet2.CPD{1}, 'cpt');
        % means = get_field(bnet2.CPD{2}, 'mean')';
        % covs = get_field(bnet2.CPD{2}, 'cov');
        % CPT3 = get_field(bnet2.CPD{3}, 'cpt');

        % figure(1)
        % clf;
        % imagesc(CPT3);

%         [LOOPERMeans] = diffusionMapTransformation(app.SavedData.BestEmission * pcaBasis, dimensionalityStream, diffusionMap, epsilons(parameterSets), ts(parameterSets), ks(parameterSets));
        LOOPERMeans = (app.SavedData.BestEmission * pcaBasis);
        LOOPERMeans(isnan(LOOPERMeans)) = 0;
        LOOPERMeans = (LOOPERMeans' - streamMean)./ streamSTD;
        LOOPERMeans = LOOPERMeans';

        allLOOPERStates = [];
        for i = 1:size(app.SavedData.BestStateMap,1)
            allLOOPERStates(i) = find(ismember(app.SavedData.BestLoopAssignments, app.SavedData.BestStateMap(i,:), 'rows'));
        end

        startDistribution = zeros(Q,1);
        startDistribution = histcounts(allLOOPERStates(trialStarts), (1:size(app.SavedData.BestModel,1)+1)-0.5);
        startDistribution = startDistribution ./ sum(startDistribution);
        % startDistribution(startState) = 1;

        trialConditions = floor((0:size(data,3)-1) / 10)+1;
%         trialConditions = mod(0:size(data,3)-1, 6)+1;

        testData = data;

        % testData = [];
        % for i = 1:size(data,3)
        %     testData(:,:,i) = pcaBasis * ((data(:,:,i) - streamMean) ./ streamSTD);
        % end

        [loglikHMM, ~, hmmStates] = mhmm_logprob(testData, CPT1, CPT3, means', covs, ones(Q,1));
        [loglikLOOPER, ~, looperStates] = mhmm_logprob(testData, startDistribution, app.SavedData.BestModel, LOOPERMeans', repmat(covs(:,:,1), [1 1 size(app.SavedData.BestModel,1)]), ones(size(app.SavedData.BestModel,1),1));

        percent = (loglikHMM - loglikLOOPER) / abs(loglikLOOPER + loglikHMM)*100

        allPercents(runNumber) = percent;

        CPT1s(:,runNumber) = CPT1;
        CPT3s(:,:,runNumber) = CPT3;
        allMeans(:,:,runNumber) = means;
    end

%     save('fmriHMM.mat', 'allPercents', 'CPT1s', 'covs', 'allMeans', 'CPT3s', 'streamSTD', 'streamMean', 'hmmStates')

    %%
    HMMConditions = [];
    allHMMStates = cell2mat(hmmStates);
    for i = 1:size(means,1)
        thisIndicies = find(allHMMStates == i);
        thisTrials = floor((thisIndicies - 1) / size(data,2)) + 1;
        thisCondition = min(trialConditions(thisTrials));

        if length(thisIndicies) == 0
            HMMConditions(i) = 1;
        else
            HMMConditions(i) = thisCondition;
        end
    end

    LOOPERConditions = [];
    % allLOOPERStates = cell2mat(looperStates);
    for i = 1:size(LOOPERMeans,1)
        thisIndicies = find(allLOOPERStates == i);
        thisTrials = floor((thisIndicies - 1) / size(data,2)) + 1;
        thisCondition = min(trialConditions(thisTrials));

        if length(thisIndicies) == 0
            LOOPERConditions(i) = 1;
        else
            LOOPERConditions(i) = thisCondition;
        end
    end

    figure(2);
    clf;
    hold on;
    for i = 1:size(testData,3)
        thisData = testData(:,:,i);

        h = plot3(thisData(1,:), thisData(2,:), thisData(3,:), 'k');
        h.Color(4) = 0.1;
    end
    scatter3(means(:,1), means(:,2), means(:,3), 'rx');


    %%

    movies = unique(eventData.movie);
    subjects = unique(eventData.subjID);

    subjects(subjects == 171734) = [];
    subjects(subjects == 103010) = [];
    subjects(subjects == 748662) = [];
    subjects(subjects == 175540) = [];

    percentCorrect = [];
    for i = 9

        moviePoolCorrect = find(contains(eventData.movie, movies{i}) & cell2mat(eventData.respAcc) == 0 & ...
                    eventData.movement == 0 & contains(eventData.matfile, 'wholeBrainPower'));
        moviePoolIncorrect = find(contains(eventData.movie, movies{i}) & cell2mat(eventData.respAcc) == 2 & ...
                    eventData.movement == 0 & contains(eventData.matfile, 'wholeBrainPower'));

        percentCorrect(i) = length(moviePoolIncorrect);
    end

    %%

    MOVIE_ID = 9;

    correctCondition = find(ismember(eventData.subjID, subjects) & contains(eventData.movie, movies{MOVIE_ID}) & cell2mat(eventData.respAcc) == 0 & ...
                    eventData.movement == 0 & contains(eventData.matfile, 'wholeBrainPower'));
    incorrectCondition = find(ismember(eventData.subjID, subjects) & contains(eventData.movie, movies{MOVIE_ID}) & cell2mat(eventData.respAcc) == 2 & ...
                    eventData.movement == 0 & contains(eventData.matfile, 'wholeBrainPower'));

    %%

    SHOULD_VALIDATE = true;
    VALDIATION_AMOUNT = 0.5;
    numBootstraps = 60;
    totalTrials = 10;
    PRE_MOVIE_FRAMES = 15;
    POST_MOVIE_FRAMES = 10;

    numBootstraps = 60;

    currentTime = clock;
    rng(floor(currentTime(6)*1000));

    conditionIDs = {};

    for i = 1:2
        if i == 1
            useIDs = correctCondition;
        else
            useIDs = incorrectCondition;
        end

        validationCount = floor(length(useIDs)*VALDIATION_AMOUNT);
        trainCount = length(useIDs) - validationCount;

        if SHOULD_VALIDATE
            conditionIDs{i} = useIDs(1:trainCount);
        else
            conditionIDs{i} = useIDs(trainCount+1:end);
        end
    end

    disp(['Building data...']);

    getDataFromIDs;

    %%

    trialData = convertToCell(allTrials);
    
    concatData = [];
    for i = 1:length(allTrials)
        concatData = [concatData allTrials{1}];
    end
    
    [~, ~, ~, ~, dataExplained] = pca(concatData);

    rawData = convertToCell(allTrials);

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

    startIndex = 20;
    
    USE_INDICIES = 1:size(rawData,1);
    
        %%

    for modelType = 0:1
        %% Reconstruct
        if modelType == 1
            bestStates = allHMMStates;
            bestModel = CPT3;
%             [~, originalData] = diffusionMapTransformation(((means' .* streamSTD) + streamMean)', dimensionalityStream, diffusionMap, epsilons(parameterSets), ts(parameterSets), ks(parameterSets));
%             clusterMeans = originalData * pcaBasis';
%             clusterMeans = ((means' .* streamSTD) + streamMean)' * pcaBasis';
        else
            bestStates = allLOOPERStates;
            bestModel = app.SavedData.BestModel;
%             clusterMeans = app.SavedData.BestEmission;
        end
        
        originalTrialLength = size(allTrials{1},2);
        trialLength = size(app.SavedData.FinalStream,1) / length(allTrials);

%         finalStream = app.SavedData.RawData(:,finalIDs)';
        finalStream = app.SavedData.FinalStream;

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
        for i = 1:length(allTrials)
            thisTrialIDs = (i-1)*trialLength + (1:trialLength);
            reconstructedStream(:, i, :) = clusterMeans(bestStates(thisTrialIDs),USE_INDICIES)';
        end

        %% Calc correlation

        trialRsquareds = [];
        for i = 1:size(reconstructedStream,2)
            thisTrialIDs = (i-1)*trialLength + (1:trialLength);
            targetData = finalStream(thisTrialIDs,USE_INDICIES)';

            thisReconstruction = squeeze(reconstructedStream(:,i, :));

            trialRsquareds(i) = corr(thisReconstruction(:), targetData(:));
        end

        meanTrialRsquared = nanmedian(trialRsquareds(:))

        if modelType == 1
            hmmDistribution = trialRsquareds;
        else
            looperDistribution =  trialRsquareds;
        end


    end

    if exist('hmmDistribution') && exist('looperDistribution')

        figure(1);
        clf;
        hold on;
        ksdensity(hmmDistribution(:));
        ksdensity(looperDistribution(:));
    
        p = ranksum(mean(hmmDistribution,2), mean(looperDistribution,2))

        looperR2 = median(looperDistribution);
        [~,looperNumPCs] = min(abs(cumsum(dataExplained) - looperR2*100));
    end


    %%
    
    for modelType = 0:1
        %% Simulate traces
        simulateTime = size(allTrials{1},2) - startIndex;

        if modelType == 1
            bestStates = allHMMStates;
            bestModel = CPT3;
%             [~, originalData] = diffusionMapTransformation(((means' .* streamSTD) + streamMean)', dimensionalityStream, diffusionMap, epsilons(parameterSets), ts(parameterSets), ks(parameterSets));
%             clusterMeans = originalData * pcaBasis';
%             clusterMeans = ((means' .* streamSTD) + streamMean)' * pcaBasis';
        else
            bestStates = allLOOPERStates;
            bestModel = app.SavedData.BestModel;
%             clusterMeans = app.SavedData.BestEmission;
        end
        
        originalTrialLength = size(allTrials{1},2);
        trialLength = size(app.SavedData.FinalStream,1) / length(allTrials);

        %%finalStream = app.SavedData.RawData(:,finalIDs)';
        finalStream = app.SavedData.FinalStream;

        clusterMeans = [];
        for i = 1:size(bestModel,1)
            thisIndicies = find(bestStates == i);

            if isempty(thisIndicies)
                clusterMeans(i,:) = zeros(1,size(finalStream,2));
            else
                clusterMeans(i,:) = mean(finalStream(thisIndicies,:), 1);
            end
        end
        
%         finalStream = app.SavedData.FinalStream;

        allClusters = [];
        for i = 1:length(trialData)
            thisStart = trialData{i}(:,startIndex);

            clusterDistances = [];
            for j = 1:size(app.SavedData.BestEmission,1)
                thisIndices = find(bestStates == j);

                stds = std(finalStream(thisIndices,:), [], 1);

                thisEmission = mean(finalStream(thisIndices,:), 1);

                clusterStream = (thisEmission - thisStart') ./ stds;           

                clusterDistances(j,:) = sum(clusterStream.^2,2);
            end

            [~, bestClusters] = min(clusterDistances, [], 1);

            allClusters(i) = bestClusters;        

        end



        reconstructedStream = [];
        for i = 1:length(trialData)            
            for j = 1:100
                currentCluster = allClusters(i);

                for t = 2:simulateTime
                    currentCluster(t) = randsample(1:size(bestModel,1), 1, true, bestModel(currentCluster(t-1),:));
                end

                reconstructedStream(:, i, :, j) = clusterMeans(currentCluster,USE_INDICIES)';
            end
        end

        %% Calc correlation

        SIM_TIME = 20;
        useTimes = startIndex:startIndex+SIM_TIME-1;

        trialRsquareds = [];
        for i = 1:size(reconstructedStream,2)
            targetData = squeeze(trialData{i}(USE_INDICIES,useTimes));

            for repeatCount = 1:size(reconstructedStream,4)
                thisReconstruction = squeeze(reconstructedStream(:,i,1:SIM_TIME, repeatCount));

                trialRsquareds(i,repeatCount) = corr(thisReconstruction(:), targetData(:));
            end
        end

        meanTrialRsquared = nanmedian(trialRsquareds(:))

        ksdensity(trialRsquareds(:));

        if modelType == 1
            hmmDistribution = trialRsquareds;
        else
            looperDistribution =  trialRsquareds;
        end


    end

    if exist('hmmDistribution') && exist('looperDistribution')
        figure(1);
        clf;
        hold on;
        ksdensity(hmmDistribution(:));
        ksdensity(looperDistribution(:));
        
        p = ranksum(mean(hmmDistribution,2), mean(looperDistribution,2))
    end
end

%% Decoding rates
decodeRateLOOPER = [];
decodeRateHMM = [];
decodeRatio = [];

for bootstrapTrial = 1:10
    for USE_HMM = 0:1
        TR = 0.72;

        movieLength = 25;
        plotLength = 42;

        plotTime = [0:plotLength-1] - (plotLength - movieLength - POST_MOVIE_FRAMES);
        plotTime = plotTime * TR;



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
                conditionOccupancies(thisState, j, trialConditions(i)) = conditionOccupancies(thisState, j, trialConditions(bootstrapID)) + 1;
            end
        end

        conditionOccupanciesStateMarginalized = squeeze(sum(conditionOccupancies,2));
        conditionOccupanciesStateMarginalized = conditionOccupanciesStateMarginalized ./ sum(conditionOccupanciesStateMarginalized,2);
        conditionOccupanciesStateMarginalized(isnan(conditionOccupanciesStateMarginalized)) = 0;


        %% Esitamte decoding rate

        loopOccupancies = conditionOccupancies ./ sum(conditionOccupancies,3);
        loopOccupancies(isnan(loopOccupancies)) = 0;

        responseTime = 37;

        F1Indices = 28:38;

        F1Indices(F1Indices > trialLength) = [];

        F1Rates = [];
        for i = 1:2
            conditionID = i;

            totalDecodeRates = [];
            for j = 1:trialLength
                loopDistribution = conditionOccupancies(:,j,i);
                loopDistribution = loopDistribution / sum(loopDistribution);

                totalDecodeRates(j) = sum(loopDistribution .* loopOccupancies(:,j,i));                
            end

            F1Rates(i) = mean(abs(totalDecodeRates(F1Indices) - 1));
        end
        F1Rate = 1 - mean(F1Rates);

        if USE_HMM == 1
            decodeRateHMM(bootstrapTrial) = F1Rate;
        else
            decodeRateLOOPER(bootstrapTrial) = F1Rate;
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
            plot(plotTime, squeeze(allTimeProbabilities(j,:,j)));
        end

        plot(plotTime(ones(2,1)*responseTime), [0 1], 'r');
        plot([0 0], [0 1], 'b');
        plot([20 20], [0 1], 'b');
    end
    
    decodeRatio(bootstrapTrial) = (decodeRateLOOPER - 0.5) / (decodeRateHMM - 0.5);
end

%%

labels = {'LOOPER Simulation R^2', 'HMM Simulation R^2', 'LOOPER Decode Rate', 'HMM Decode Rate'};
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