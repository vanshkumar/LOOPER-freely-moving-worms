
%% Load validation data FOR RNN

currentTime = clock;
rng(floor(currentTime(6)*1000));

clear time;


load('goodRNN.mat');
load('goodRNNResult');

app.SavedData = saveData;

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

    % From RNN

    NUM_TRIALS = 10;

    times = 1:220;

    useClasses = [0,1,2,3,4,5];

    trialData = permute(dynamics, [3, 1, 2]);

    allData = reshape(trialData, size(trialData,1),[]);

    [pcaBasis,~,~,~,dataExplained] = pca(allData', 'NumComponents', 20);

%     pcaBasis = eye(size(trialData,1));

%     figure(1);
%     clf;
%     colors = lines(length(useClasses));
%     hold on;

    trialCounts = [];
    finalTrials = [];
    rawTrials = [];
    trialID = 1;
    for i = 1:length(useClasses)
        thisIndices = find(classes == useClasses(i) & missTrials == 0);

%         F1ID = mod(useClasses(i), 3) + 1;
%         F2ID = floor(useClasses(i) / 3) + 1;


        for j = (1:NUM_TRIALS)+USE_VALIDATION*NUM_TRIALS
            thisData = (squeeze(trialData(:,times,thisIndices(j)))')';
            for k = 1:size(thisData,1)
                finalTrials(k,:, trialID) = decimate(thisData(k,:),2);
            end

            thisData = (squeeze(trialData(:,times,thisIndices(j)))' * pcaBasis)';
            for k = 1:size(thisData,1)
                rawTrials(k,:, trialID) = decimate(thisData(k,:),2);
            end

            trialID = trialID + 1;
        end

%         plot(mean(finalInputs(1,:,end-9:end),3), 'Color', colors(i,:));
    end

    numTrial = max(app.SavedData.TrialData);
    originalLength = size(saveData.RawData,2)/numTrial;
    trialLength = size(saveData.FinalStream,1) / numTrial;
    
    rawData = convertToCell(rawTrials(:,3:end,:));
    lastEnd = 0;
    rawTrialData = [];
    for i = 1:length(rawData)
        rawTrialData(lastEnd + (1:size(rawData{i}, 2))) = i;
    
        lastEnd = lastEnd + size(rawData{i}, 2);
    end
    rawData = mergeData(rawData);

    trialData = convertToCell(finalTrials(:,1+originalLength-trialLength:originalLength,:)) ;
    
    dataPCABAsis = pcaBasis;
catch
end

%% Load LOOPER data

load('BayesoptResultsRNN.mat')

pcCount = BayesoptResults.XAtMinObjective.pcCount;
epsilon = BayesoptResults.XAtMinObjective.epsilon;
tValue = BayesoptResults.XAtMinObjective.t;
kValue = BayesoptResults.XAtMinObjective.k;
stateCounts = BayesoptResults.XAtMinObjective.stateCounts;


decimateAmount = 2;

times = 1:900;

app.SavedData = saveData;

numTrial = max(app.SavedData.TrialData);
trialLength = size(app.SavedData.FinalStream,1) / numTrial;
plotTimes = 1:500;

if max(plotTimes) > trialLength
    plotTimes = plotTimes(1):trialLength;
end

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
%     trialConditions = mod(0:size(data,3)-1, 6)+1;

    testData = data;

    % testData = [];
    % for i = 1:size(data,3)
    %     testData(:,:,i) = (data(:,:,i) - streamMean) ./ streamSTD;
    % end

    [loglikHMM, ~, hmmStates] = mhmm_logprob(testData, CPT1, CPT3, means', covs, ones(Q,1));
    
    looperCOVs = repmat(eye(size(covs,1)) * 0.1, [1 1 size(app.SavedData.BestModel,1)]);
    
    [loglikLOOPER, ~, looperStates] = mhmm_logprob(testData, startDistribution, app.SavedData.BestModel, LOOPERMeans', looperCOVs, ones(size(app.SavedData.BestModel,1),1));

    percent = (loglikHMM - loglikLOOPER) / abs(loglikLOOPER + loglikHMM)*100

    allPercents(runNumber) = percent;
    
    CPT1s(:,runNumber) = CPT1;
    CPT3s(:,:,runNumber) = CPT3;
    allMeans(:,:,runNumber) = means;
end

save('rnnHMM.mat', 'allPercents', 'CPT1s', 'covs', 'allMeans', 'CPT3s', 'streamSTD', 'streamMean')

%%

meanPercent = mean(allPercents);
sePercent = std(allPercents) * 1.98 / sqrt(length(allPercents));
disp([num2str(meanPercent) ' +/- ' num2str(sePercent)]);

%%

HMMConditions = [];
HMMConditionWeights = [];
allHMMStates = cell2mat(hmmStates);
for i = 1:size(means,1)
    thisIndicies = find(allHMMStates == i);
    thisTrials = floor((thisIndicies - 1) / size(data,2)) + 1;
    thisConditions = histcounts(trialConditions(thisTrials), (0:6)+0.5);
    [thsiWeights, thisConditions] = sort(thisConditions, 'descend');
    
    if length(thisIndicies) == 0
        HMMConditions(i,:) = 1:2;
        HMMConditionWeights(i,:) = [1, 1];
    else
        HMMConditions(i,:) = thisConditions(1:2);
        HMMConditionWeights(i,:) = thsiWeights(1:2);
    end
end

trialF1s = [];
trialF2s = [];
for i = 1:6
    trialF1s(i) = mod(i - 1,3)+1;
    trialF2s(i) = floor((i-1)/3)+1;
end

LOOPERConditions = [];
LOOPERConditionWeights = [];
% allLOOPERStates = cell2mat(looperStates);
for i = 1:size(LOOPERMeans,1)
    thisIndicies = find(allLOOPERStates == i);
    thisTrials = floor((thisIndicies - 1) / size(data,2)) + 1;
    thisConditions = histcounts(trialConditions(thisTrials), (0:6)+0.5);
    [thsiWeights, thisConditions] = sort(thisConditions, 'descend');
    
    if length(thisIndicies) == 0
        LOOPERConditions(i,:) = 1:2;
        LOOPERConditionWeights(i,:) = [1, 1];
    else
        LOOPERConditions(i,:) = thisConditions(1:2);
        LOOPERConditionWeights(i,:) = thsiWeights(1:2);
    end
end

%% tSNE

hmmGraph = digraph(CPT3, 'omitselfloops');
looperGraph = digraph(app.SavedData.BestModel, 'omitselfloops');


POWER = 8;
EXAGGERATION = 2;
PERPLEXITY = 40;
NUM_DIMENSIONS = 3;

HMMDist = @(i,j) precomputedDistanceFunction(i,j,mk_stochastic(CPT3 + CPT3')^POWER);
LOOPERDist = @(i,j) precomputedDistanceFunction(i,j,mk_stochastic(app.SavedData.BestModel + app.SavedData.BestModel')^POWER);

options = statset('TolFun',eps);

indicies = repmat((1:size(CPT3,1))', [1 100]);
hmmY = tsne(indicies, 'Algorithm', 'exact', 'Distance', HMMDist, 'NumDimensions', NUM_DIMENSIONS, 'Perplexity', PERPLEXITY, 'Options', options, 'Exaggeration', EXAGGERATION);
indicies = repmat((1:size(app.SavedData.BestModel,1))', [1 100]);
looperY = tsne(indicies, 'Algorithm', 'exact', 'Distance', LOOPERDist, 'NumDimensions', NUM_DIMENSIONS, 'Perplexity', PERPLEXITY, 'Options', options, 'Exaggeration', EXAGGERATION);

clear 'lines';
lineColors = [1 0.6 0; 1 0 0.6; 0 1 0.6; 0.6 1 0; 0.6 0 1; 0 0.6 1];
lineColors = 0.1 + lineColors*0.9;

lineThetas = deg2rad([40; 160; 280; 320; 80; 200] + 10);
lineS = 0.8;
lineV = 0.8;

figure(2);
clf
subplot(1,2,1);
if size(hmmY,2) == 2
%     scatter(hmmY(:,1), hmmY(:,2));
    plot(hmmGraph, 'xData', hmmY(:,1), 'yData', hmmY(:,2))
else
    HMMAlphas = HMMConditionWeights(:,1) ./ sum(HMMConditionWeights,2);
    HMMColors = [];
    HMMColors(:,:,1) = lineColors(HMMConditions(:,1),:) .* repmat(HMMAlphas, [1 3]);
    HMMColors(:,:,2) = lineColors(HMMConditions(:,2),:) .* (repmat(1 - HMMAlphas, [1 3]));
    
    HMMColors = HMMColors(:,:,1) + HMMColors(:,:,2);
    
    HMMThetas = [];
    HMMThetas(:,1) = exp(1i*lineThetas(HMMConditions(:,1),:)) .* HMMAlphas;
    HMMThetas(:,2) = exp(1i*lineThetas(HMMConditions(:,2),:)) .* (1 - HMMAlphas);
    
    allAngles = angle(HMMThetas(:,1) + HMMThetas(:,2));
    allAngles = wrapTo2Pi(allAngles) / (2*pi);
    
    HMMColors = squeeze(hsv2rgb([allAngles, repmat(lineS, [size(allAngles,1),1]), repmat(lineV, [size(allAngles,1),1])]));
    
%     scatter3(hmmY(:,1), hmmY(:,2), hmmY(:,3));
    scatter3(hmmY(:,1), hmmY(:,2), hmmY(:,3), 72, HMMColors, 'filled');
    hold on;
    for i = 1:size(hmmGraph.Edges,1)        
        startNode = hmmGraph.Edges.EndNodes(i,1);
        
        fromPoint = hmmY(startNode, :);
        toPoint = hmmY(hmmGraph.Edges.EndNodes(i,2), :);
        
%         points = [(fromPoint + toPoint)/2; toPoint];
        points = [fromPoint; (fromPoint + toPoint)/2];
        
        h = plot3(points(:,1), points(:,2), points(:,3), 'LineWidth', 2, 'Color', HMMColors(startNode,:));
        h.Color(4) = 0.3;
    end
    
    [~, startNode] = max(CPT1);
    scatter3(hmmY(startNode,1), hmmY(startNode,2), hmmY(startNode,3), 150, 'rx');
    
    title('HMM');
%     plot(hmmGraph, 'xData', hmmY(:,1), 'yData', hmmY(:,2), 'zData', hmmY(:,3))
end
subplot(1,2,2);

if size(looperY,2) == 2
%     scatter(looperY(:,1), looperY(:,2));
%     plot(looperGraph, 'xData', looperY(:,1), 'yData', looperY(:,2))
else
    LOOPERAlphas = LOOPERConditionWeights(:,1) ./ sum(LOOPERConditionWeights,2);
    LOOPERColors = [];
    LOOPERColors(:,:,1) = lineColors(LOOPERConditions(:,1),:) .* repmat(LOOPERAlphas, [1 3]);
    LOOPERColors(:,:,2) = lineColors(LOOPERConditions(:,2),:) .* (repmat(1 - LOOPERAlphas, [1 3]));
    
    LOOPERColors = LOOPERColors(:,:,1) + LOOPERColors(:,:,2);
    
    LOOPERThetas = [];
    LOOPERThetas(:,1) = exp(1i*lineThetas(LOOPERConditions(:,1),:)) .* LOOPERAlphas;
    LOOPERThetas(:,2) = exp(1i*lineThetas(LOOPERConditions(:,2),:)) .* (1 - LOOPERAlphas);
    
    allAngles = angle(LOOPERThetas(:,1) + LOOPERThetas(:,2));
    allAngles = wrapTo2Pi(allAngles) / (2*pi);
    
    LOOPERColors = squeeze(hsv2rgb([allAngles, repmat(lineS, [size(LOOPERThetas,1),1]), repmat(lineV, [size(LOOPERThetas,1),1])]));
    
    
    scatter3(looperY(:,1), looperY(:,2), looperY(:,3), 72, LOOPERColors, 'filled');
    hold on;
    for i = 1:size(looperGraph.Edges,1)        
        startNode = looperGraph.Edges.EndNodes(i,1);
        
        fromPoint = looperY(startNode, :);
        toPoint = looperY(looperGraph.Edges.EndNodes(i,2), :);
        
%         points = [(fromPoint + toPoint)/2; toPoint];
        points = [fromPoint; (fromPoint + toPoint)/2];
        
        h = plot3(points(:,1), points(:,2), points(:,3), 'LineWidth', 2, 'Color', LOOPERColors(startNode,:));
        h.Color(4) = 0.3;
    end
    title('LOOPER');
%     plot(looperGraph, 'xData', looperY(:,1), 'yData', looperY(:,2), 'zData', looperY(:,3))
end
% plot(hmmGraph, 'xData', Y(:,1), 'yData', Y(:,2))

%% Reconstruction

for USE_HMM = 0:1
    if USE_HMM
        bestStates = allHMMStates;
        bestModel = CPT3;
    else
        bestStates = allLOOPERStates;
        bestModel = app.SavedData.BestModel;
    end

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
    [~,looperNumPCs] = min(abs(cumsum(dataExplained) - looperR2*100));
end




%% Simulate traces

for USE_HMM = 0:1
    startIndex = 40;
    
    simulateTime = size(finalTrials,2) - startIndex;

    if USE_HMM
        bestStates = allHMMStates;
        bestModel = CPT3;
        % clusterMeans = ((means' .* streamSTD) + streamMean)' * pcaBasis';
    else
        bestStates = allLOOPERStates;
        bestModel = app.SavedData.BestModel;
        % clusterMeans = app.SavedData.BestEmission;
    end

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

    [tempData, procesedTrialData, procesedTrialSwitches] = preprocessData(rawData, [], 0, [], 0, rawTrialData, 0, true, saveData.PreprocessData.Smoothing, saveData.PreprocessData.ZScore, saveData.PreprocessData.DelayTime, saveData.PreprocessData.DelayCount, saveData.DataMean, saveData.DataSTD);
    finalStream = app.SavedData.FinalStream;

    allClusters = [];
    for i = 1:length(trialData)
        thisStart = procesedTrialData{i}(:,startIndex);

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

    SIM_TIME = 40;
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
                conditionOccupancies(thisState, j, trialConditions(i)) = conditionOccupancies(thisState, j, trialConditions(bootstrapID)) + 1;
            end
        end

        conditionOccupanciesStateMarginalized = squeeze(sum(conditionOccupancies,2));
        conditionOccupanciesStateMarginalized = conditionOccupanciesStateMarginalized ./ sum(conditionOccupanciesStateMarginalized,2);
        conditionOccupanciesStateMarginalized(isnan(conditionOccupanciesStateMarginalized)) = 0;


        %% Esitamte decoding rate

        loopOccupancies = conditionOccupancies ./ sum(conditionOccupancies,3);
        loopOccupancies(isnan(loopOccupancies)) = 0;

        F1Start = 10;
        F2Start = 35;

        F1Indices = 20:30;
        F2Indices = 45:55;

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
        figure(10+USE_HMM);
        clf;
        hold on;
        for j = 1:length(unique(trialConditions))
            plot(squeeze(allTimeProbabilities(j,:,j)));
        end

        plot(ones(2,1)*F1Start, [0 1], 'r');
        plot(ones(2,1)*F2Start, [0 1], 'r');
    end
    
    decodeRatio(bootstrapTrial) = (decodeRateLOOPER - 0.25) / (decodeRateHMM - 0.25);
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
