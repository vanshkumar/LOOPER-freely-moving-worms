USE_BROKEN_RNN = max(classes) < 9;

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

NUM_TRIALS = 10;

times = 1:220;

if ~USE_BROKEN_RNN
    useClasses = [1,2,3,6,7,8];
else
    useClasses = [0,1,2,3,4,5];
end

trialData = permute(dynamics, [3, 1, 2]);

allData = reshape(trialData, size(trialData,1),[]);

[pcaBasis,~] = pca(allData', 'NumComponents', 20);

% pcaBasis = eye(size(size(trialData,1)));

% figure(1);
% clf;
% colors = lines(length(useClasses));
% hold on;

finalTrials = [];
finalInputs = [];
finalOutputs = [];
finalTargets = [];
trialCounter = 1;
rawTrials = [];
for i = 1:length(useClasses)
    thisIndices = find(classes == useClasses(i) & missTrials == 0);
    
    if IS_VALIDATION
        thisIndices(1:NUM_TRIALS) = [];
    end
    
    for j = 1:NUM_TRIALS
        thisData = (squeeze(trialData(:,times,thisIndices(j)))' * pcaBasis)';
        for k = 1:size(thisData,1)
            finalTrials(k,:,trialCounter) = decimate(thisData(k,:),2);
        end
        finalInputs(:,:,trialCounter) = (inputs(:,times(1):2:times(end),thisIndices(j)));
        finalOutputs(1,:,trialCounter) = (double(outputs(thisIndices(j),times(1):2:times(end))));
        finalTargets(:,:,trialCounter) = (targets(:,times(1):2:times(end),thisIndices(j)));
        rawTrials(:,:,trialCounter) = thisData;
        
        trialCounter = trialCounter + 1;
    end
    
%     plot(mean(finalInputs(1,:,end-9:end),3), 'Color', colors(i,:));
end


%%

% First load RNN or Monkey LOOPER data

PLOT_START = 45;
PLOT_END = 179;
F1Start = 50;
F2Start = 110;

STATE_SMOOTH = 2;
FLUX_CUTOFF = 3;
MINIMUM_STATE_TIME = 0 *(2+1);
IS_RNN = 1;
decimateTime = 2;
SHOULD_VALIDATE = IS_VALIDATION;

app.SavedData = saveData;



numTrial = max(app.SavedData.TrialData);
originalLength = size(saveData.RawData,2)/numTrial;
trialLength = size(app.SavedData.FinalStream,1) / numTrial;
times = 1:900;

if max(times) > trialLength
    times = times(1):trialLength;
end

trialIndicies = repmat(times, [1, numTrial]);
trialIndicies = trialIndicies + kron((0:numTrial-1)*trialLength, ones(size(times)));

finalStream = app.SavedData.FinalStream;

[trialPCABasis, ~] = pca(finalStream(trialIndicies,:), 'NumComponents',3);
% trialPCABAsis = pcaBasis;

loopStarts = (0:numTrial-1)*trialLength+1;

finalStream(loopStarts,:) = nan;


if SHOULD_VALIDATE
%     matSize = size(allBootstrappedFiringRates);
%     matData = permute(allBootstrappedFiringRates, [1, 4, 2, 3, 5]);
%     allTrials = reshape(matData, [matSize(1), matSize(4), matSize(2) * matSize(3) * matSize(5)]);
% 
%     inputData = permute(allInputs, [1, 4, 2, 3, 5]);
%     allInputs = reshape(inputData, [1, matSize(4), matSize(2) * matSize(3) * matSize(5)]);
% 
% 
%     finalTrials = [];
%     finalInputs = [];
%     for i = 1:size(allTrials,1)
%         for j = 1:size(allTrials,3)
%             finalTrials(i,:,j) = decimate(squeeze(allTrials(i,:,j)),decimateTime);
% 
%             if i == 1
%                 finalInputs(1,:,j) = decimate(squeeze(allInputs(1,:,j)),decimateTime);
%             end
%         end
%     end
    
    rawData = convertToCell(finalTrials);

    lastEnd = 0;
    rawTrialData = [];
    for i = 1:length(rawData)
        rawTrialData(lastEnd + (1:size(rawData{i}, 2))) = i;

        lastEnd = lastEnd + size(rawData{i}, 2);
    end
    
    rawData = mergeData(rawData);

%     rawData = app.SavedData.RawData;
%     rawTrialData = app.SavedData.TrialData;

    [tempData, trialData, procesedTrialSwitches] = preprocessData(rawData, [], 0, [], 0, rawTrialData, 0, true, saveData.PreprocessData.Smoothing, saveData.PreprocessData.ZScore, saveData.PreprocessData.DelayTime, saveData.PreprocessData.DelayCount, saveData.DataMean, saveData.DataSTD);
    badIndicies = [procesedTrialSwitches procesedTrialSwitches+1];
    badIndicies(badIndicies > size(tempData,2)) = [];
    tempData(:,badIndicies) = [];
    
%     tempData = app.SavedData.FinalStream';

%     figure(4);
%     clf;
%     plot(tempData(1,:));
%     hold on;
%     plot(app.SavedData.FinalStream(:,1)');

    
    clusterDistances = [];
    finalStream = app.SavedData.FinalStream;
    for i = 1:size(app.SavedData.BestEmission,1)
        thisLoopPosition = app.SavedData.BestLoopAssignments(i,:);
        thisIndices = find(ismember(app.SavedData.BestStateMap, thisLoopPosition, 'rows'));
        
        stds = std(finalStream(thisIndices,:), [], 1);
        
        thisEmission = mean(finalStream(thisIndices,:), 1);
        
        clusterStream = (thisEmission - tempData') ./ stds;           
        
        clusterDistances(i,:) = sum(clusterStream.^2,2);% ./ sqrt(length(thisIndices));
    end
    
%     bestClusters = [];
%     clusterDistances = pdist2(app.SavedData.BestEmission, tempData');
    [~, bestClusters] = min(clusterDistances, [], 1);

    loopIDs = app.SavedData.BestLoopAssignments(bestClusters,1);
    
    MINIMUM_STATE_TIME = 5;
else    
    loopIDs = app.SavedData.BestStateMap(:,1);
end
if MINIMUM_STATE_TIME > 0
    loopIDs = colfilt(loopIDs, [MINIMUM_STATE_TIME 1], 'sliding', @mode);
end

if IS_RNN
    if USE_BROKEN_RNN
        loopOrder = [5, 8, 7, 6, 1, 2, 4, 3];
    else
        loopOrder = [6 5 1 2 4 3 7 8];
    end

    lineStyles = {'b', 'g', 'r', 'b.', 'g.', 'r.'};
    lineColors = {'b', 'g', 'r', 'b', 'g', 'r'};
    
    originalTime = 1:size(dynamics,1);
    decimateAmount = 2;
else
    loopOrder = [6, 5, 4, 1, 2, 3, nan];
%     loopOrder = [4, 6, 5, 3, 1, 2, nan];

    lineStyles = {'r', 'g', 'b', 'r.', 'g.', 'b.'};
    lineColors = {'r', 'g', 'b', 'r', 'g', 'b'};
    
    decimateAmount = decimateTime;
    originalTime = [-2500 10000];
    originalTime = originalTime(1)/1000 : 0.01 : originalTime(2)/1000;
end

plotTime = originalTime((1:decimateAmount:trialLength*decimateAmount) + (originalLength - trialLength)*decimateAmount);

trialLoopIDs = [];
trialConditionIDs = [];

figureHandle = figure(figureNumber);
figureHandle.Renderer='Painters';
clf;
hold on;
legendLines = [];
for i = 1:numTrial
    trialIDs = (0:trialLength-1)+(i-1)*trialLength+1;
    
    if IS_RNN
        conditionID = floor((i-1)/10)+1;
    else
        conditionID = mod((i-1),6)+1;
    end
    
    IDs = loopIDs(trialIDs);
    IDs(IDs == 0) = length(loopOrder);
    
    h = plot(plotTime, loopOrder(IDs) + normrnd(0,0.03,size(trialIDs)), lineStyles{conditionID});
    h.Color(4) = 0.2;
    legendLines(conditionID) = h;
    
    trialLoopIDs(:,i) = loopOrder(IDs);
    trialConditionIDs(i) = conditionID;
    
    transitionCheck = colfilt(loopOrder(IDs), [1 2], 'sliding', @mean);
    transitionCheck(transitionCheck == loopOrder(IDs)) = 0;
    transitionPoints = find(transitionCheck);
    transitionPoints = [transitionPoints, transitionPoints+1];
    transitionPoints(transitionPoints > trialLength) = [];
    
    transitionData = nan(size(loopOrder(IDs)));
    transitionData(transitionPoints) = loopOrder(IDs(transitionPoints));
    
    h = plot(plotTime, transitionData + normrnd(0,0.03,size(trialIDs)), lineColors{conditionID});
    h.Color(4) = 0.2;
end
yLimit = ylim;
plot(ones(1,2)*F1Start, yLimit, 'k');
plot(ones(1,2)*F2Start, yLimit, 'k');
xlim([PLOT_START PLOT_END]);

% legend(legendLines, {'10Hz -> 6Hz', '18Hz -> 10Oz', '34Hz -> 26Hz', '10Hz -> 18Hz', '18Hz -> 26Hz', '34HZ -> 44Hz'});
% legend(legendLines, {'10Hz -> 5Hz', '20Hz -> 15Hz', '40Hz -> 30Hz', '10Hz -> 15Hz', '20Hz -> 25Hz', '40HZ -> 50Hz'});
