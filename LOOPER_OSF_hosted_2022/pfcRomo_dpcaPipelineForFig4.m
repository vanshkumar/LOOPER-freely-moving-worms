
% Broken RNN
load('F:\Dropbox\ConservationOfAgentDynamics\WorkingMemoryTask\Experiments\workingMemoryCustomBroken\lstm\1\networkTrace.mat');

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
sum(missTrials) / length(missTrials)

% From RNN

NUM_TRIALS = 20;

times = 1:180;

useClasses = [0,1,2,3,4,5];

trialData = permute(dynamics, [3, 1, 2]);

allData = reshape(trialData, size(trialData,1),[]);

% [pcaBasis,~] = pca(allData', 'NumComponents', 50);

pcaBasis = eye(size(trialData,1));

% figure(1);
% clf;
% colors = lines(length(useClasses));
% hold on;

trialCounts = [];
finalTrials = [];
for i = 1:length(useClasses)
    thisIndices = find(classes == useClasses(i) & missTrials == 0);
    
    F1ID = mod(useClasses(i), 3) + 1;
    F2ID = floor(useClasses(i) / 3) + 1;
    
    
    for j = 1:NUM_TRIALS
        thisData = (squeeze(trialData(:,times,thisIndices(j)))' * pcaBasis)';
        for k = 1:size(thisData,1)
            trialCounts(k,F1ID, F2ID) = NUM_TRIALS;
            finalTrials(k,F1ID, F2ID, :, j) = decimate(thisData(k,:),2);
        end
    end
    
%     plot(mean(finalInputs(1,:,end-9:end),3), 'Color', colors(i,:));
end

%%
% clear all
% load data_romo_eLife.mat  % this file is produced by pfcRomo_preprocess.m

% firingRates array is stored in the file in compressed sparse format.
% This line is de-compressing it.
% firingRates = reshape(full(firingRates_sparse), firingRates_size);
% spikeRasters = reshape(full(spikeRasters_sparse), spikeRasters_size);

% Neuron selection criteria used in the eLife paper
% D = size(trialNum,1);
% minN = min(reshape(trialNum(:,:,:), D, []), [], 2);
% meanFiringRate = mean(reshape(firingRatesAverage, D, []), 2);
% n = find(minN >= 5 & meanFiringRate < 50);
% 
% firingRates = firingRates(n,:,:,:,:);
% firingRatesAverage = firingRatesAverage(n,:,:,:);
% spikeRasters = spikeRasters(n,:,:,:,:);
% trialNum = trialNum(n,:,:);

firingRates = finalTrials;
firingRatesAverage = mean(finalTrials,5);
trialNum = trialCounts;

% IMPORTANT NOTE: This yields 788 neurons, instead of 832 as reported in
% the eLife paper. This discrepancy is because we had a mistake in the
% preprocessing script and selected as neurons some auxilliary channels
% that are actually not neurons. This does not influence the results in any
% substantial way, because these auxilliary channels are mostly silent. To
% obtain the same number of units as in the paper, run
% pfcRomo_preprocess.m with electrodeNum = 8 instead of electrodeNum = 7.
% This will yield 832 units.

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Stimulus', 'Decision', 'Condition-independent', 'Interaction'};
decodingClasses = {[1 1; 2 2; 3 3; 4 4; 5 5; 6 6], [1 2; 1 2; 1 2; 1 2; 1 2; 1 2], [], [1 2; 3 4; 5 6; 7 8; 9 10; 11 12]};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

%% Cross-validation to find lambda

% This takes some time (around 4*10 min on my laptop) and produces 
% optimalLambda = 2.5629e-06;

optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
    'combinedParams', combinedParams, ...
    'numComps', [10 10 10 10], ...
    'numRep', 10, ...
    'filename', 'tmp_optimalLambdasFig4.mat');

%% dPCA (with regularization and noise cov)

Cnoise = dpca_getNoiseCovariance(firingRatesAverage, firingRates, trialNum, ...
    'type', 'averaged');
  
[W,V,whichMarg] = dpca(firingRatesAverage, 50, ...
    'combinedParams', combinedParams, 'lambda', optimalLambda, 'Cnoise', Cnoise);

%%

doZScore = 0;
numDelays = 0;
delayTime = 0;

X = firingRatesAverage(:,:)';
Xcen = bsxfun(@minus, X, mean(X));
Z = Xcen * W;
dataDim = size(firingRatesAverage);

XTrial = firingRates(:,:)';
XTrialscen = bsxfun(@minus, XTrial, mean(X));
ZTrial = XTrialscen * W;
dataDimTrial = size(firingRates);

componentsToPlot = 1:50;

Zfull = reshape(Z(:,componentsToPlot)', [length(componentsToPlot) dataDim(2:end)]);
ZTrialfull = reshape(ZTrial(:,componentsToPlot)', [length(componentsToPlot) dataDimTrial(2:end)]);

colors = jet(7);
colors = colors([1 2 4 5 6 7],:);
colors(3,:) = colors(3,:)*0.6;
colors(4,:) = [238 210 2]/256;

colors = colors([1, 3, 6], :);

% colors = {'r', 'g', 'b'};

stimMarg = find(strcmp(margNames, 'Stimulus'));
decisionMarg = find(strcmp(margNames, 'Decision'));
interactionMarg = find(strcmp(margNames, 'Interaction'));

topStim = find(whichMarg == stimMarg, 1);
topDecision = find(whichMarg == decisionMarg, 1);
topInteraction = find(whichMarg == interactionMarg, 1);

plotComponents = [topStim, topDecision];
times = 1:90;

F1Time = 50/2;
F2Time = 110/2;

allTraces = [];
for f=1:size(Zfull,2)
    for d = 1:size(Zfull,3)
        thisTrace = squeeze(Zfull(plotComponents, f, d, times));        
        if doZScore
            thisTrace = zscore(thisTrace, [], 2);
        end
        thisTrace = delayEmbed(thisTrace, numDelays, delayTime);
        
        allTraces = [allTraces thisTrace];
    end
end

if size(allTraces,1) < 3
    pcaBasis = eye(size(allTraces,1),3);
else
    [pcaBasis, ~] = pca(allTraces', 'NumComponents', 6);
end

figureHandle = figure(3);
figureHandle.Renderer='Painters';
clf;
hold on;
for f=[1,2,3]
    for d = 1:size(Zfull,3)
        if d == 1
            lineType = '.';
        else
            lineType = '-';
        end
        
        thisTrace = squeeze(Zfull(plotComponents, f, d, times));
        
%         [pcaBasis, ~] = pca(thisTrace', 'NumComponents', 3);

%         plot3(thisTrace'*pcaBasis(:,1), thisTrace'*pcaBasis(:,2), thisTrace'*pcaBasis(:,3), lineType, 'color', colors(f,:), 'LineWidth', 2)
        plot(thisTrace'*pcaBasis(:,1), lineType, 'color', colors(f,:), 'LineWidth', 2)
        
        for j = 1:size(ZTrialfull,5)
            thisTrace = squeeze(ZTrialfull(plotComponents, f, d, times, j));
            
            h = plot(thisTrace'*pcaBasis(:,1), lineType, 'color', colors(f,:));
            h.Color(4) = 0.3;
        end
%         scatter3(thisTrace'*pcaBasis(:,1), thisTrace'*pcaBasis(:,2), thisTrace'*pcaBasis(:,3), 32, 1:size(thisTrace,2))
    end
end

% Stats
traces = squeeze(ZTrialfull(plotComponents(1), 2:3, 1:2, times, :));

probabilities = [];
isDifferent = [];
for i = 1:length(times)
    distribution1 = squeeze(traces(1,1,i,:));
    distribution2 = squeeze(traces(2,1,i,:));
    distribution3 = squeeze(traces(1,2,i,:));
    distribution4 = squeeze(traces(2,2,i,:));
    
    allDistributions = [[distribution1; distribution3], [distribution2; distribution4]];
    
    probabilities(i) = kruskalwallis(allDistributions, [], 'off');
    
%     [probabilities(i), isDifferent(i)] = ranksum(allDistributions(:,1), allDistributions(:,2));
end

isDifferent = double(probabilities < 0.01);

isDifferent(isDifferent == 0) = nan;

yLimit = ylim;

plot(isDifferent*(yLimit(1) + 1), 'k', 'LineWidth',3);

plot([F1Time F1Time], yLimit);
plot([F2Time F2Time], yLimit);



%% Stop here for minimum processing
% 
% explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
%      'combinedParams', combinedParams, ...
%      'Cnoise', Cnoise, 'numOfTrials', trialNum);
% 
% dpca_plot(firingRatesAverage, W, V, @dpca_plot_romo, ...
%     'whichMarg', whichMarg,                 ...
%     'time', time,                           ...
%     'timeEvents', timeEvents,               ...
%     'timeMarginalization', 3,               ...
%     'ylims', [150 150 400 150],             ...
%     'legendSubplot', 16,                    ...
%     'marginalizationNames', margNames,      ...
%     'explainedVar', explVar,                ...
%     'marginalizationColours', margColours);
% 
% %% decoding part
% 
% % with 100 iterations this takes around 10*100/60 = 17 min on my laptop
% 
% accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
%     'lambda', optimalLambda, ...
%     'combinedParams', combinedParams, ...
%     'decodingClasses', decodingClasses, ...
%     'noiseCovType', 'averaged', ...
%     'numRep', 5, ...        % increase to 100
%     'filename', 'tmp_classification_accuracy.mat');
% 
% dpca_classificationPlot(accuracy, [], [], [], decodingClasses)
% 
% % with 100 iterations and 100 shuffles this takes 100 times longer than the
% % above function, i.e. 17*100/60 = 28 hours (on my laptop). Be careful.
% 
% accuracyShuffle = dpca_classificationShuffled(firingRates, trialNum, ...
%     'lambda', optimalLambda, ...
%     'combinedParams', combinedParams, ...
%     'decodingClasses', decodingClasses, ...
%     'noiseCovType', 'averaged', ...
%     'numRep', 5, ...        % increase to 100
%     'numShuffles', 20, ...  % increase to 100 (takes a lot of time)
%     'filename', 'tmp_classification_accuracy.mat');
% 
% dpca_classificationPlot(accuracy, [], accuracyShuffle, [], decodingClasses)
% 
% componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);
% 
% dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
%     'explainedVar', explVar, ...
%     'marginalizationNames', margNames, ...
%     'marginalizationColours', margColours, ...
%     'whichMarg', whichMarg,                 ...
%     'time', time,                        ...
%     'timeEvents', timeEvents,               ...
%     'timeMarginalization', 3,           ...
%     'legendSubplot', 16,                ...
%     'componentsSignif', componentsSignif);
