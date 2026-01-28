function [correlation, constraints] = ComapreToMarkovFunction(pcCount, epsilon, tValue, kValue, stateCounts, ranges)

   try   

    if exist('ranges') && ~isempty(ranges)
        pcCount = round(pcCount * (ranges.pcCountMax - ranges.pcCountMin)) + ranges.pcCountMin;
        epsilon = (epsilon * (ranges.epsilonMax - ranges.epsilonMin)) + ranges.epsilonMin;
        tValue = round(tValue * (ranges.tMax - ranges.tMin)) + ranges.tMin;
        kValue = round(kValue * (ranges.kMax - ranges.kMin)) + ranges.kMin; 
        stateCounts = round(stateCounts * (ranges.stateCountsMax - ranges.stateCountsMin)) + ranges.stateCountsMin;
    end
    
        %% Load validation data FOR MONKEY
    
    currentTime = clock;
    rng(floor(currentTime(6)*1000));

    clear time;

    romoAnalysisPipeline_forFig3
%     load('data_romo_eLife.mat');
%     firingRates = full(firingRates_sparse);
%     firingRates = reshape(firingRates, firingRates_size);
    forValidation = 1;
    monkeySetupForFig3
    % load('monkeyDataForFigure3')
    
    load('monkeyDataForFigure3.mat');

    app.SavedData = saveData;

    numTrial = max(saveData.TrialData);
    originalLength = size(saveData.RawData,2)/numTrial;
    trialLength = size(saveData.FinalStream,1) / numTrial;
    
    rawTrials = [];
    for i = 1:size(finalTrials,3)
        rawTrials(:,:,i) = dataPCABAsis * finalTrials(:,:,i);
    end

    trialData = convertToCell(rawTrials(:,1+originalLength-trialLength:originalLength,:)) ;

%     rawData = convertToCell(finalTrials);
% 
%     lastEnd = 0;
%     rawTrialData = [];
%     for i = 1:length(rawData)
%         rawTrialData(lastEnd + (1:size(rawData{i}, 2))) = i;
% 
%         lastEnd = lastEnd + size(rawData{i}, 2);
%     end
% 
%     rawData = mergeData(rawData);
% 
%     [tempData, trialData, procesedTrialSwitches] = preprocessData(rawData, [], 0, [], 0, rawTrialData, 0, true, saveData.PreprocessData.Smoothing, saveData.PreprocessData.ZScore, saveData.PreprocessData.DelayTime, saveData.PreprocessData.DelayCount, saveData.DataMean, saveData.DataSTD);
            
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

    dimensionalityStream = app.SavedData.FinalStream;

    pcaBasis = eye(size(dimensionalityStream,2));

    [pcaBasis, ~, ~, ~, explained] = pca(dimensionalityStream);
    pcaCutoff = pcCount;
    pcaBasis = pcaBasis(:, 1:pcaCutoff);

    dimensionalityStream = dimensionalityStream * pcaBasis;

    rawStream = dimensionalityStream;

    [rawStream, diffusionMap] = diffusionMapReduction(dimensionalityStream, epsilon, tValue, kValue);


    streamMean = mean(rawStream', 2);
    streamSTD = std(rawStream', [], 2);

    data = [];
    for i = 1:numTrial
        finalindices = trialStarts(i):trialStarts(i)+trialLength-1;

%         data(:,:,i) = [ones(size(rawStream,2), 1)*(1+rand(1)*0.1) (rawStream(finalindices,1:size(rawStream,2))' - streamMean) ./ streamSTD ones(size(rawStream,2), 1)*(1+rand(1)*0.1)];
        data(:,:,i) = (rawStream(finalindices,1:size(rawStream,2))' - streamMean) ./ streamSTD;
    end

    O = size(data,1);          %Number of coefficients in a vector 
    T = size(data,2);         %Number of vectors in a sequence 
    nex = size(data,3);        %Number of sequences 
    M = 1;          %Number of mixtures 
    Q = stateCounts;          %Number of states 
    % Q = 110;          %Number of states 
    cov_type = 'diag';

    %% Build HMM
    
    disp(['Trainning HMM with params: pcs=' num2str(pcCount) ' eps=' num2str(epsilon) ' t=' num2str(tValue) ' k=' num2str(kValue) ' #=' num2str(stateCounts)] );
    

    allPercents = [];
    CPT1s = [];
    CPT3s = [];
    allMeans = [];
    for runNumber = 1:1

        matrixPrior = zeros(Q,Q);
        % matrixPrior = diag(ones(Q-1,1),1);
        % matrixPrior(Q,1) = 1;
        % matrixPrior = kron(eye(3),matrixPrior) * 5000/Q/3;
        % matrixPrior = matrixPrior + eye(size(matrixPrior)) * 5000/Q/6;


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
            covariates(:,:,i) = eye(O)*0.1;
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
        bnet.CPD{3} = tabular_CPD(bnet, 3, 'CPT', transmat0);
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
        trueMeans = (get_field(bnet2.CPD{2}, 'mean') .* streamSTD + streamMean)';
        covs = get_field(bnet2.CPD{2}, 'cov');
        CPT3 = get_field(bnet2.CPD{3}, 'cpt');

        figure(1);
        clf;
        subplot(2,1,1);
        imagesc(means);
        subplot(2,1,2);
        imagesc(app.SavedData.BestEmission);


        %% Compare models
        % CPT1 = get_field(bnet2.CPD{1}, 'cpt');
        % means = get_field(bnet2.CPD{2}, 'mean')';
        % covs = get_field(bnet2.CPD{2}, 'cov');
        % CPT3 = get_field(bnet2.CPD{3}, 'cpt');

        % figure(1)
        % clf;
        % imagesc(CPT3);

        LOOPERMeans = app.SavedData.BestEmission * pcaBasis;
        LOOPERMeans(isnan(LOOPERMeans)) = 0;
        LOOPERMeans = (LOOPERMeans' - streamMean)./ streamSTD;
        LOOPERMeans = LOOPERMeans';

        allLOOPERStates = [];
        for i = 1:size(app.SavedData.BestStateMap,1)
            allLOOPERStates(i) = find(ismember(app.SavedData.BestLoopAssignments, app.SavedData.BestStateMap(i,:), 'rows'));
        end

        startDistribution = zeros(size(app.SavedData.BestModel,1),1);
        startDistribution = histcounts(allLOOPERStates(trialStarts), (1:size(app.SavedData.BestModel,1)+1)-0.5);
        startDistribution = startDistribution ./ sum(startDistribution);
        % startDistribution(startState) = 1;

        trialConditions = floor((0:size(data,3)-1) / 10)+1;
        trialConditions = mod(0:size(data,3)-1, 6)+1;

        testData = data;

        % testData = [];
        % for i = 1:size(data,3)
        %     testData(:,:,i) = (data(:,:,i) - streamMean) ./ streamSTD;
        % end

        [loglikHMM, ~, hmmStates] = mhmm_logprob(testData, CPT1, CPT3, means', covs, ones(Q,1));
        allHMMStates = cell2mat(hmmStates);

        looperCOVs = repmat(eye(size(covs,1)) * 0.1, [1 1 size(app.SavedData.BestModel,1)]);

        [loglikLOOPER, ~, looperStates] = mhmm_logprob(testData, startDistribution, app.SavedData.BestModel, LOOPERMeans', looperCOVs, ones(size(app.SavedData.BestModel,1),1));

        percent = (loglikHMM - loglikLOOPER) / abs(loglikLOOPER + loglikHMM)*100

        allPercents(runNumber) = percent;

        CPT1s(:,runNumber) = CPT1;
        CPT3s(:,:,runNumber) = CPT3;
        allMeans(:,:,runNumber) = means;
    end


    %% Simulate traces

    simulateTime = size(finalTrials,2) - startIndex;

    bestStates = allHMMStates;
    bestModel = CPT3;

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


    meanTrialRsquared = nanmean(trialRsquareds(:))

%     figureHandle = figure(9);
%     figureHandle.Renderer='Painters';
%     clf;
%     ksdensity(trialRsquareds(:));
    
    
    correlation = 1 - abs(meanTrialRsquared);
    constraints = -1;
   catch ME
    warning(ME.message);
    disp([ME.stack(1).name ' Line: ' num2str(ME.stack(1).line)]);
    correlation = 1;
    constraints = 1;
   end
end

