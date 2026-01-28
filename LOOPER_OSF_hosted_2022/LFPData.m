%% Loops

load('\\192.168.1.206\LabData\adeeti\ConnorCollab\GL13_archExp\awake_GL13_LFP.mat');

awakeData = meanSubData;
awakeData(info.noiseChannels,:,:) = [];

load('LFPAwake2LoopsMoreTrialsFinal');

% load('\\192.168.1.206\LabData\adeeti\ConnorCollab\GL13_archExp\highIso_GL13_LFP_DiskStation_Aug-27-1112-2020_CaseConflict.mat')

% highIsoData = meanSubData;
% highIsoData(info.noiseChannels,:,:) = [];

% allData = {awakeData, highIsoData};
allData = {awakeData};

%%

concatenatedData = [];
for j = 1:length(allData)
    for i = 1:size(allData{j},1)
        concatenatedData = [concatenatedData squeeze(allData{j}(:,i,:))];
    end
end

%%

[pcaBasis,~,~,~,explained] = pca(concatenatedData');
percentExplained = cumsum(explained);
numDimensions = find(percentExplained > 95, 1);
% numDimensions = length(percentExplained);
% pcaBasis = eye(numDimensions);

percentExplained = percentExplained(numDimensions)

%%

SHOULD_VALIDATE = true;
totalTrials = 300; %50 for model fits
numBootstraps = 5;
validationPercent = 0.5;
NUM_TRIALS = size(allData{1},2);
USE_TIME = 500:2000;
DECIMATE_AMOUNT = 10;

rng(1);

useIDs = randsample(1:NUM_TRIALS, NUM_TRIALS);

currentTime = clock;
% rng(floor(currentTime(6)*1000));

rng(15); %Make sure we see a variety of trial types

allTrials = [];
trialIndex = 1;
for i = 1:length(allData)
    validIDs = 1:NUM_TRIALS;
    
    allIDs = 1:floor(length(validIDs)*validationPercent);
    validationIDs = floor(length(validIDs)*validationPercent)+1:length(useIDs);
    
    for j = 1:totalTrials
        if ~SHOULD_VALIDATE
            bootstrapIDs = randsample(length(allIDs), numBootstraps, true);
        else
            bootstrapIDs = randsample(length(validationIDs), numBootstraps, true);
        end
        thisIDs = useIDs(validIDs(bootstrapIDs));
        
        thisData = zeros(size(allData{i},1), length(decimate(USE_TIME, DECIMATE_AMOUNT)));
        for k = 1:length(thisIDs)
            for l = 1:size(allData{i},1)
                thisData(l,:) = thisData(l,:) + decimate(squeeze(allData{i}(l,thisIDs(k), USE_TIME)), DECIMATE_AMOUNT)';
            end
        end
        thisData = thisData / length(thisIDs);

        allTrials(:,:,trialIndex) = (thisData' * pcaBasis(:,1:numDimensions))';

        trialIndex = trialIndex + 1;
    end
end

%%
SHOULD_VALIDATE = true;
MINIMUM_STATE_TIME = 0 *(2+1);

CHECK_FRAMES = 85:129;
RINGING_LOOP = 5;

CHECK_FRAMES_INITIAL = 45:80;
INITIAL_LOOP = 2;

app.SavedData = saveData;

numTrial = max(saveData.TrialData);
reducedSize = size(saveData.TrialData,2) - size(saveData.FinalStream,1);
reducedSize = reducedSize / numTrial;
trialLength = size(saveData.FinalStream,1) / numTrial;
% times = 1:trialLength;

% if max(times) > trialLength
%     times = times(1):trialLength;
% end

% trialIndicies = repmat(times, [1, numTrial]);
% trialIndicies = trialIndicies + kron((0:numTrial-1)*trialLength, ones(size(times)));

finalStream = saveData.FinalStream;

if size(saveData.ClusterMeans,2) > 2
    [trialPCABasis, ~] = pca(saveData.ClusterMeans, 'NumComponents',3);
else
    trialPCABasis = eye(size(saveData.ClusterMeans,2));
end

loopStarts = [0, saveData.RawTrialSwitches] - [0:numTrial-1]*reducedSize + 1;

longestTrial = max(diff([loopStarts size(finalStream,2)]));

finalStream(loopStarts,:) = nan;

loopStarts(end+1) = size(finalStream,1);

if SHOULD_VALIDATE    
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

    
    clusterDistances = [];
    finalStream = saveData.FinalStream;
    for i = 1:size(saveData.BestEmission,1)
        thisLoopPosition = saveData.BestLoopAssignments(i,:);
        thisIndices = find(ismember(saveData.BestStateMap, thisLoopPosition, 'rows'));
        
        stds = std(finalStream(thisIndices,:), [], 1);
        
        thisEmission = mean(finalStream(thisIndices,:), 1);
        
        clusterStream = (thisEmission - tempData') ./ stds;           
        
        clusterDistances(i,:) = sum(clusterStream.^2,2);
    end
    
    [~, bestClusters] = min(clusterDistances, [], 1);

    loopIDs = saveData.BestLoopAssignments(bestClusters,1);
    
    MINIMUM_STATE_TIME = 10;
else    
    loopIDs = saveData.BestStateMap(:,1);
end
if MINIMUM_STATE_TIME > 0
    loopIDs = colfilt(loopIDs, [MINIMUM_STATE_TIME 1], 'sliding', @mode);
end

loopOrder = [1, 3, 1, 2, 5, 6, 7, 8, 9, 10, 11, 12];
% loopOrder = [4, 3, 5, 1, 2, 6, 7, 8, 9, 10, 11, 12];

clear lines;
lineColors = lines(length(loopOrder));

% trialLoopIDs = [];
% trialConditionIDs = [];

ringingness = [];
intialness = [];

figureHandle = figure(1);
figureHandle.Renderer='Painters';
clf;
hold on;
legendLines = [];
for i = 1:length(trialData)    
    conditionID = 1;%floor((i-1)/10)+1;
%     trialIDs = loopStarts(i):loopStarts(i+1);
    trialIDs = (1:longestTrial) + (i-1)*longestTrial;
    
    IDs = loopIDs(trialIDs);
    IDs(IDs == 0) = 1;
    
    plotTime = 1:length(IDs);
    
    h = plot(plotTime, loopOrder(IDs) + normrnd(0,0.03,size(trialIDs)), 'Color', lineColors(conditionID,:));
    h.Color(4) = 0.1;
    legendLines(conditionID) = h;
    
    ringingness(i) = sum(IDs(CHECK_FRAMES) == RINGING_LOOP);
    intialness(i) = sum(IDs(CHECK_FRAMES_INITIAL) == INITIAL_LOOP);
    
%     trialLoopIDs(:,i) = loopOrder(IDs);
%     trialConditionIDs(i) = conditionID;
    
%     transitionCheck = colfilt(loopOrder(IDs), [1 2], 'sliding', @mean);
%     transitionCheck(transitionCheck == loopOrder(IDs)) = 0;
%     transitionPoints = find(transitionCheck);
%     transitionPoints = [transitionPoints, transitionPoints+1];
%     transitionPoints(transitionPoints > trialLength) = [];
%     
%     transitionData = nan(size(loopOrder(IDs)));
%     transitionData(transitionPoints) = loopOrder(IDs(transitionPoints));
    
%     h = plot(plotTime, transitionData + normrnd(0,0.03,size(trialIDs)), lineColors{conditionID});
%     h.Color(4) = 0.2;
end

% legend(legendLines, {'10Hz -> 6Hz', '18Hz -> 10Oz', '34Hz -> 26Hz', '10Hz -> 18Hz', '18Hz -> 26Hz', '34HZ -> 44Hz'});
% legend(legendLines, {'10Hz -> 5Hz', '20Hz -> 15Hz', '40Hz -> 30Hz', '10Hz -> 15Hz', '20Hz -> 25Hz', '40HZ -> 50Hz'});
    

lowInitial = find(intialness < 10);
highRinging = find(ringingness > length(CHECK_FRAMES)*0.9);

PLOT_COUNT = 10;
CAP_COUNTS = 2;
useLow = lowInitial(randsample(length(lowInitial), CAP_COUNTS, 'false'));
useHigh = highRinging(randsample(length(highRinging), CAP_COUNTS, 'false'));

useMiddle = round(linspace(CAP_COUNTS, length(ringingness) - CAP_COUNTS, PLOT_COUNT - 2*CAP_COUNTS +2));
useMiddle = useMiddle(2:end-1);

[~, plotOrder] = sort(ringingness);

offsetValue = 1;

figureHandle = figure(2);
figureHandle.Renderer='Painters';
clf;
hold on;
legendLines = [];
% for i = 1:PLOT_COUNT
for i = 1:length(ringingness)
%     if i <= CAP_COUNTS
%         thisTrial = useLow(i);
%     elseif i > PLOT_COUNT - CAP_COUNTS
%         thisTrial = useHigh(i - PLOT_COUNT + CAP_COUNTS);
%     else
%         thisTrial = plotOrder(useMiddle(i - CAP_COUNTS));
%     end
    
    thisTrial = plotOrder(i);
    
    
    conditionID = 1;%floor((i-1)/10)+1;
%     trialIDs = loopStarts(i):loopStarts(i+1);
    trialIDs = (1:longestTrial) + (thisTrial-1)*longestTrial;
    
    thisData = tempData(1, trialIDs);
    
    IDs = loopIDs(trialIDs);
    IDs(IDs == 0) = 1;
    
    plotTime = 1:length(IDs);
    
    plotColor(plotTime, thisData + (i-1) * offsetValue, loopOrder(IDs), lines(8));
end

USE_TRIALS = [5 7 26 103 49 56 60 74 113 174 226 270 289 293 295];

figureHandle = figure(3);
figureHandle.Renderer='Painters';
clf;
hold on;
legendLines = [];
for i = 1:length(USE_TRIALS)
    thisTrial = plotOrder(USE_TRIALS(i));
    
    
    conditionID = 1;%floor((i-1)/10)+1;
%     trialIDs = loopStarts(i):loopStarts(i+1);
    trialIDs = (1:longestTrial) + (thisTrial-1)*longestTrial;
    
    thisData = tempData(1, trialIDs);
    
    IDs = loopIDs(trialIDs);
    IDs(IDs == 0) = 1;
    
    plotTime = 1:length(IDs);
    
    plotColor(plotTime, thisData + (i-1) * offsetValue, loopOrder(IDs), lines(8));
end

%%

plotReconstruction;

numTrial = 1;

figure(1);
clf;
hold on;
trialIDs = loopStarts(numTrial):loopStarts(numTrial+1);
plot(processedStream(1,trialIDs));
plot(reconstructedStream(1,trialIDs));

%% Run Sims

SHOULD_VALIDATE = true;
totalTrials = 50;
numBootstraps = 5;
validationPercent = 0.5;
NUM_TRIALS = size(allData{1},2);
USE_TIME = 500:2000;
DECIMATE_AMOUNT = 10;

rng(1);

useIDs = randsample(1:NUM_TRIALS, NUM_TRIALS);

currentTime = clock;
rng(floor(currentTime(6)*1000));

allTrials = [];
rawTrials = [];
trialIndex = 1;
for i = 1:length(allData)
    validIDs = 1:NUM_TRIALS;
    
    allIDs = 1:floor(length(validIDs)*validationPercent);
    validationIDs = floor(length(validIDs)*validationPercent)+1:length(useIDs);
    
    for j = 1:totalTrials
        if ~SHOULD_VALIDATE
            bootstrapIDs = randsample(length(allIDs), numBootstraps, true);
        else
            bootstrapIDs = randsample(length(validationIDs), numBootstraps, true);
        end
        thisIDs = useIDs(validIDs(bootstrapIDs));
        
        thisData = zeros(size(allData{i},1), length(decimate(USE_TIME, DECIMATE_AMOUNT)));
        for k = 1:length(thisIDs)
            for l = 1:size(allData{i},1)
                thisData(l,:) = thisData(l,:) + decimate(squeeze(allData{i}(l,thisIDs(k), USE_TIME)), DECIMATE_AMOUNT)';
            end
        end
        thisData = thisData / length(thisIDs);

        allTrials(:,:,trialIndex) = (thisData' * pcaBasis(:,1:numDimensions))';
        rawTrials(:,:,trialIndex) = thisData';

        trialIndex = trialIndex + 1;
    end
end

rawStream = saveData.FinalStream(:,1:size(saveData.RawData,1));

clusterMeans = [];
clusterSTD = [];
clusterLengths = [];
meanTimes = [];
stateLookup = [];
for i = 1:size(saveData.BestLoopAssignments,1)
    loopID = saveData.BestLoopAssignments(i,1);
    phaseID = saveData.BestLoopAssignments(i,2);
    
    thisIndices = find(saveData.BestStateMap(:,1) == loopID & saveData.BestStateMap(:,2) == phaseID);
    
    trialTimes = mod(thisIndices, trialLength)+1;
    meanTimes(i,:) = mode(trialTimes);
    clusterLengths(i,:) = length(thisIndices);

    clusterMeans(i,:) = nanmean(rawStream(thisIndices, :), 1);
    clusterSTD(i,:) = nanstd(rawStream(thisIndices, :), [], 1);
    
    stateLookup(loopID, phaseID) = i;
end

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


clusterDistances = [];
finalStream = saveData.FinalStream;
for i = 1:size(saveData.BestEmission,1)
    thisLoopPosition = saveData.BestLoopAssignments(i,:);
    thisIndices = find(ismember(saveData.BestStateMap, thisLoopPosition, 'rows'));

    stds = std(finalStream(thisIndices,:), [], 1);
    thisEmission = mean(finalStream(thisIndices,:), 1);

    clusterStream = (thisEmission - tempData') ./ stds;           

    clusterDistances(i,:) = sum(clusterStream.^2,2);
end

[~, bestClusters] = min(clusterDistances, [], 1);



trialTime = floor(length(bestClusters)/size(allTrials,3));

START_FRAME = 63;

bestDataClusters = [];
for i = 1:size(saveData.BestEmission,1)
    thisIDs = find(ismember(saveData.BestStateMap, saveData.BestLoopAssignments(i,:), 'row'));
    
    bestDataClusters(thisIDs) = i;
end

bestTransitions = [];
for i = 1:size(saveData.BestModel,1)
    thisIndices = find(bestDataClusters == i);
    thisIndices(thisIndices == length(bestDataClusters)) = [];
    
    for j = 1:size(saveData.BestModel,2)
        bestTransitions(i,j) = sum((bestDataClusters(thisIndices+1) == j));
    end
end
bestTransitions = bestTransitions ./ sum(bestTransitions,2);

NUM_SIMS = 100;

reconstructedTrials = [];
for i = 1:size(allTrials,3)
    startFrame = START_FRAME + 1 - (size(allTrials,2) - trialTime) + (i-1)*trialTime;
    startState = bestClusters(startFrame);
    
    for k = 1:NUM_SIMS
        stateTrial = nan(START_FRAME-1,1);
        currentState = startState;
        for t = START_FRAME:trialTime
            stateTrial(t) = currentState;

%             currentState = randsample(1:saveData.BestStateCount, 1, true, bestTransitions(currentState,:));
            currentState = randsample(1:saveData.BestStateCount, 1, true, saveData.BestModel(currentState,:));
        end
        
        reconstructedTrials(i,1:START_FRAME,:,k) = nan(START_FRAME, size(clusterMeans,2));
        reconstructedTrials(i,START_FRAME:trialTime,:,k) = clusterMeans(stateTrial(START_FRAME:trialTime),:);
    end
end

meanReconstructions = nanmean(reconstructedTrials, 4);
stdReconstructions = nanstd(reconstructedTrials, [], 4);

%%

goodChannels = setdiff(1:64, info.noiseChannels);

USE_TIME = 50;

PLOT_TRIAL = 37;

PLOT_CHANNELS = [52, 25, 49, 20];

projectedData = squeeze(allTrials(:,1:START_FRAME+USE_TIME,PLOT_TRIAL))';
projectedData = projectedData * pcaBasis(:,1:numDimensions)';

rawData = squeeze(rawTrials(1:START_FRAME+USE_TIME,:,PLOT_TRIAL));

allRawData = squeeze(rawTrials(1:START_FRAME+USE_TIME,:,:));
meanRawData = mean(allRawData,3);
stdRawData = std(allRawData,[],3);

reconstructionMeans = squeeze(meanReconstructions(PLOT_TRIAL,START_FRAME:START_FRAME+USE_TIME,:)) * pcaBasis(:,1:numDimensions)';
reconstructionSTDs = squeeze(stdReconstructions(PLOT_TRIAL,START_FRAME:START_FRAME+USE_TIME,:)) * pcaBasis(:,1:numDimensions)';

figureHandle = figure(1);
figureHandle.Renderer='Painters';
clf;

plotChannels = [];
for i = 1:length(PLOT_CHANNELS)
    plotChannel = find(goodChannels == PLOT_CHANNELS(i));
    subplot(length(PLOT_CHANNELS),1,i);
    hold on;
%     plot(rawData(:,plotChannel), 'k', 'LineWidth', 1);
    plot(projectedData(:,plotChannel), 'k', 'LineWidth', 1);
%     plotShadedCI(meanRawData(:,plotChannel), stdRawData(:,plotChannel), 1:START_FRAME+USE_TIME, 'k');
    plotShadedCI(reconstructionMeans(:,plotChannel), reconstructionSTDs(:,plotChannel), START_FRAME:START_FRAME+USE_TIME, 'r');
    title(['Channel ' num2str(PLOT_CHANNELS(i))]);
    
    plotChannels(i) = plotChannel;
end

%%

reconstructionPCABasis = pcaBasis(:,1:numDimensions)';
reconstructionPlotIndices = 1:1000;

plotReconstruction;

% 
% %% setting up outline
% outline = imread(['F:\\Dropbox\\MethodsPaperCode\\Data\\', 'MouseBrainAreas.png']);
% darknessOutline = 1;
% mmInGridX = 2.75;
% mmInGridY = 5;
% PIXEL_TO_MM = (2254 - 503)/2;
% BREGMA_PIXEL_X = 5520;
% BREGMA_PIXEL_Y = 3147;
% bregmaOffsetPixelX = info.bregmaOffsetX * PIXEL_TO_MM;
% bregmaOffsetPixelY = info.bregmaOffsetY * PIXEL_TO_MM;
% overlayWindow = [[0 round(mmInGridX*PIXEL_TO_MM)]; [0 round(mmInGridY*PIXEL_TO_MM)]];
% overlayWindow(1,:) = overlayWindow(1,:) - round(mmInGridX*PIXEL_TO_MM) + BREGMA_PIXEL_X - bregmaOffsetPixelX;
% overlayWindow(2,:) = overlayWindow(2,:) + BREGMA_PIXEL_Y + bregmaOffsetPixelY;
% outline = imgaussfilt(outline, 4);
% alpha = outline(:,:,1) / max(max(outline(:,:,1))) * darknessOutline;
% %%%%%%
% %%%%%%
% 
% figureHandle = figure(2);
% figureHandle.Renderer='Painters';
% clf;
% hold on;
% 
% imagesc(alpha);
% 
% ha = axes;
% % plot([1 size(info.gridIndicies,2) size(info.gridIndicies,2) 1 1], [1 1 size(info.gridIndicies,1) size(info.gridIndicies,1) 1], 'k');
% for i = 1:length(PLOT_CHANNELS)
%     gridIndex = find(info.gridIndicies == PLOT_CHANNELS(i));
%     [y, x] = ind2sub(size(info.gridIndicies),gridIndex);
%     y = size(info.gridIndicies,1)-y+1;
%     
% %     [y, x] = ([y, x] - 1)*PIXEL_TO_MM;
%     
%     plot(x,size(info.gridIndicies,1)-y+1,'rx');
% end
% xlim([0 size(info.gridIndicies,2)+1]);
% ylim([0 size(info.gridIndicies,1)+1]);
% 
% 
% a1 = axes;
% hold on;
% a1.Position = ha.Position;
% h = imshow(outline);
% set(h, 'AlphaData', alpha);
% xlim(overlayWindow(1,:));
% ylim(overlayWindow(2,:));
% overlayAspectRatio = (overlayWindow(1,2) - overlayWindow(1,1))/(overlayWindow(2,2) - overlayWindow(2,1));
% dataAspectRatio = ha.PlotBoxAspectRatio(1) / ha.PlotBoxAspectRatio(2);
% a1.DataAspectRatioMode = 'manual';
% a1.DataAspectRatio = [overlayAspectRatio/dataAspectRatio, 1, 1];
% 
% %%



%%
% for i = 1:size(reconstructedTrials,3)
%     h = plot(squeeze(reconstructedTrials(PLOT_TRIAL,:,1,i)), 'r');
%     h.Color(4) = 0.3;
% end



allCorrelations = [];
for i = 1:size(allTrials,3)
    thisTrial = squeeze(allTrials(:,START_FRAME:START_FRAME+USE_TIME,i));
    thisTrial = thisTrial' * pcaBasis(:,1:numDimensions)';
    
    for k = 1:NUM_SIMS
        thisSim = squeeze(reconstructedTrials(i,START_FRAME:START_FRAME+USE_TIME,:,k))';
        thisSim = thisSim' * pcaBasis(:,1:numDimensions)';
        
        allCorrelations(i,k) = corr(thisTrial(:), thisSim(:));
    end
end

[~,bestTrials] = sort(sum(allCorrelations,2), 'descend');

figureHandle = figure(3);
figureHandle.Renderer='Painters';
clf;
hold on;
ksdensity(allCorrelations(:));


%%

rawTraces = [];
for i = 1:50
	rawTraces(i,:,:) = saveData.RawData(1:3,(1:151)+(i-1)*151);
end

rawMeans = squeeze(mean(rawTraces));
rawSTDs = squeeze(std(rawTraces));

totalTrials = 50;
rawTrialLength = size(saveData.RawData,2)/totalTrials;
finalTrialLength = size(saveData.FinalStream,1)/totalTrials;

finalStreamOffset = rawTrialLength - finalTrialLength;

clf;
hold on;

plotShadedCI(rawMeans(1, finalStreamOffset:end), rawSTDs(1, finalStreamOffset:end));
plotShadedCI(rawMeans(2, finalStreamOffset:end), rawSTDs(1, finalStreamOffset:end), [], 'r');
plotShadedCI(rawMeans(3, finalStreamOffset:end), rawSTDs(1, finalStreamOffset:end), [], 'b');
% plot(rawMeans(finalStreamOffset:end));
yAxes = ylim;
plot([1 1]*(50-finalStreamOffset), yAxes);

% complexForm = hilbert(saveData.RawData(1,:));
% 
% plot(angle(complexForm(1:end)));

%%



clear 'reconstructionPCABasis';
plotReconstruction;

backProjection = pcaBasis(:,1:numDimensions)';

plotTimes = [67 91]+0;
% plotTimes = [72 96];

figureHandle = figure(3);
figureHandle.Renderer='Painters';
clf;
hold on;

goodChannels = setdiff(1:64, info.noiseChannels);

useCLim = [-0.05 0.1];

for i = 1:length(plotTimes)
    thisTime = plotTimes(i);
    
    thisStates = allIDs((1:finalTrialLength:length(allIDs)) + thisTime -1);
%     thisState = allIDs(thisTime);
%     thisState = mode(thisStates);
    
%     saveData.BestLoopAssignments(thisState,:)
    
    thisEmissions = saveData.BestEmission(thisStates,1:numDimensions) * backProjection;
    
    thisEmission = mean(thisEmissions);
    
    thisGrid = nan(size(info.gridIndicies'));
    for j = 1:length(thisEmission)
        thisChannel = goodChannels(j);
        
        
        gridIndex = find(info.gridIndicies == thisChannel);
        [y, x] = ind2sub(size(info.gridIndicies),gridIndex);
        
        thisGrid(x,y) = thisEmission(j);
    end
    
    subplot(1,length(plotTimes),i);
    
    
    imagesc(thisGrid');
    caxis(useCLim);
    colorbar;
end



bestPhase = 0.05;
derivativeSign = -1;

% complexForm = hilbert(saveData.RawData(1,:));
% checkValues = angle(complexForm(1:end));
% [deltaPeaks,bestIndices] = findpeaks(-abs(angleDiff(checkValues, bestPhase))); %Use for phases

checkValues = saveData.RawData(1,:);
[deltaPeaks,bestIndices] = findpeaks(-abs(checkValues - bestPhase));

bestIndices(abs(deltaPeaks) > 0.5) = [];

badIndices = [];
for i = 1:length(bestIndices)
    plusSign = sign(checkValues(bestIndices(i)+1) - checkValues(bestIndices(i)));
    minusSign = sign(checkValues(bestIndices(i)-1) - checkValues(bestIndices(i)));
    
    if plusSign == minusSign
        badIndices(end+1) = i;
    end
end
bestIndices(badIndices) = [];

if derivativeSign ~= 0
    badIndices = [];
    for i = 1:length(bestIndices)
        thisSign = sign(checkValues(bestIndices(i)+1) - checkValues(bestIndices(i)-1));

        if thisSign ~= derivativeSign
            badIndices(end+1) = i;
        end
    end
    bestIndices(badIndices) = [];
end

clf;
hold on;
plot(checkValues);
scatter(bestIndices, checkValues(bestIndices));

trialStarts = (1:rawTrialLength:size(saveData.RawData,2));
initialRange = (42:75) + finalStreamOffset;
ringingRange = (82:101) + finalStreamOffset;

initialIndices = kron(trialStarts, ones(1,length(initialRange))) + repmat(initialRange, [1 totalTrials]);
ringingIndices = kron(trialStarts, ones(1,length(ringingRange))) + repmat(ringingRange, [1 totalTrials]);


figureHandle = figure(3);
figureHandle.Renderer='Painters';
clf;
hold on;
plot(checkValues);

figureHandle = figure(4);
figureHandle.Renderer='Painters';
clf;
hold on;

goodChannels = setdiff(1:64, info.noiseChannels);

useCLim = [-0.1 0.1];

thisEmissionMeans = [];
thisEmissionSTDs = [];
pc1Means = [];
pc1STDs = [];
for i = 1:length(plotTimes)
    thisTime = plotTimes(i);
    
    if i == 1
        thisStates = find(ismember(1:size(saveData.RawData,2), initialIndices) & ismember(1:size(saveData.RawData,2), bestIndices));
    else
        thisStates = find(ismember(1:size(saveData.RawData,2), ringingIndices) & ismember(1:size(saveData.RawData,2), bestIndices));
    end
    
    thisEmissions = saveData.RawData(:,thisStates)' * backProjection;
    
    thisEmission = mean(thisEmissions);
    
    thisEmissionMeans(i,:) = thisEmission;
    thisEmissionSTDs(i,:) = std(thisEmissions);
    
    pc1Means(i) = mean(checkValues(thisStates));
    pc1STDs(i) = std(checkValues(thisStates));
    
    thisGrid = nan(size(info.gridIndicies'));
    for j = 1:length(thisEmission)
        thisChannel = goodChannels(j);
        
        
        gridIndex = find(info.gridIndicies == thisChannel);
        [y, x] = ind2sub(size(info.gridIndicies),gridIndex);
        
        thisGrid(x,y) = thisEmission(j);
    end
    
    figure(3);
    
    if i == 1
        scatter(thisStates, checkValues(thisStates), 'r');
    else
        scatter(thisStates, checkValues(thisStates), 'b');
    end
    


    figure(4);
    
    subplot(1,length(plotTimes),i);
    
    
    imagesc(thisGrid');
    caxis(useCLim);
    colorbar;
end


dPrimes = (thisEmissionMeans(1,:) - thisEmissionMeans(2,:)) ./ sqrt((thisEmissionSTDs(1,:).^2 + thisEmissionSTDs(2,:).^2)/2);

pcDPrime = (pc1Means(1) - pc1Means(2)) ./ sqrt((pc1STDs(1).^2 + pc1STDs(2).^2)/2)

figureHandle = figure(5);
figureHandle.Renderer='Painters';
clf;

thisGrid = nan(size(info.gridIndicies'));
for j = 1:length(thisEmission)
    thisChannel = goodChannels(j);


    gridIndex = find(info.gridIndicies == thisChannel);
    [y, x] = ind2sub(size(info.gridIndicies),gridIndex);

    thisGrid(x,y) = dPrimes(j);
end

imagesc(thisGrid');
colorbar;
