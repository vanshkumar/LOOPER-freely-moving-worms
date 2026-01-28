
phaseCounts = 100;
phases = (0:phaseCounts-1) / phaseCounts * 2 * pi;
sigmas = 5 * exp(-((phases-pi).^2)/(2*(pi/3)^2));
sigmas = sigmas + 5;

sizeX = 100;
sizeY = 100;

energyLandscape = zeros(sizeY, sizeX);

radius = min(sizeX, sizeY) / 2 - 20;

theta = 0;

positionX = [50];
positionY = [50];
depth = [.5];

thetaSign = 1;
for j = 1:1
    for i = 1:length(phases)   
        theta = phases(i);

        currentX = mod(round(positionX(j)+cos(thetaSign*theta)*radius-1),sizeX)+1;
        currentY = mod(round(positionY(j)+sin(thetaSign*theta)*radius-1),sizeY)+1;

        x = (1:sizeX) - currentX;
        y = (1:sizeY) - currentY;
        [X, Y] = meshgrid(x, y);

        sigma = sigmas(i);
        kernel = depth(j) * exp(-(X.^2 + Y.^2)/(2*sigma^2));

        energyLandscape = bsxfun(@max, energyLandscape, kernel);
    end
end

energyLandscape = energyLandscape / max(energyLandscape(:));

figure(1);
clf;
imagesc(energyLandscape);

%%

startSigma = 1;
endSigma = 1;

maxPoints = 600;

% energyLandscape = drawValley(458, 59, 800, 800, startSigma, endSigma, maxPoints);
% energyLandscape = energyLandscape .* eraseValley(900, 800, 100, 800);
% 
% tempValley = drawValley(900, 800, 100, 800, startSigma, endSigma, maxPoints);
% energyLandscape = bsxfun(@max, energyLandscape, tempValley .* eraseValley(100, 800, 500, 100));
% 
% tempValley = drawValley(50, 885, 500, 150, startSigma, endSigma, maxPoints);
% energyLandscape = bsxfun(@max, energyLandscape, tempValley .* eraseValley(458, 59, 800, 800));
% 
% energyLandscape = 1 - energyLandscape;
% 
% figure(1);
% clf;
% imagesc(energyLandscape);

sizeX = 100;
sizeY = 100;



% radius = min(sizeX, sizeY) / 2 - 10;

theta = 0;

positionX = [50, 38];
positionY = [50, 50];
depth = [1, .7];
radius = [35, 20];
sigmas = [8, 4];

wellEnergy = zeros(sizeY, sizeX, size(depth,2));

thetaSign = 1;
%waitHandle = parfor_progressbar(maxPoints * 2, 'Building energy landscape');
for j = 1:size(depth,2)
    for i = 0:maxPoints    
        theta = i / maxPoints * (2 * pi * 1.1);

        currentX = mod(round(positionX(j)+cos(thetaSign*theta)*radius(j)-1),sizeX)+1;
        currentY = mod(round(positionY(j)+sin(thetaSign*theta)*radius(j)-1),sizeY)+1;

        x = (1:sizeX) - currentX;
        y = (1:sizeY) - currentY;
        [X, Y] = meshgrid(x, y);

        sigma = startSigma + (i/maxPoints)^0.5 * (endSigma - startSigma);
        kernel = depth(j) * exp(-(X.^2 + Y.^2)./(2*sigmas(j)^2));

        wellEnergy(:,:,j) = bsxfun(@max, wellEnergy(:,:,j), kernel);

        %waitHandle.iterate(1);
    end
end
%close(waitHandle);
% energyLandscape = 1 - (1 - wellEnergy(:,:,1)) .* (1 - wellEnergy(:,:,2));
energyLandscape = sum(wellEnergy,3);
energyLandscape = energyLandscape / max(energyLandscape(:));

figure(1);
clf;
colormap(jet)
imagesc(energyLandscape)
% contour(energyLandscape);


%%



interpolateEnergyLandscape = @(x,y) interp2(1:sizeX,1:sizeY,energyLandscape,x,y);

discreteSize = 0.01;

gradUx = @(x,y) (interpolateEnergyLandscape(x - discreteSize/2,y) - interpolateEnergyLandscape(x + discreteSize/2,y)) / discreteSize;
gradUy = @(x,y) (interpolateEnergyLandscape(x,y - discreteSize/2) - interpolateEnergyLandscape(x,y + discreteSize/2)) / discreteSize;

x1 = (2:1:sizeX-1);
y1 = (2:1:sizeY-1);
[X, Y] = meshgrid(x1, y1);

gaussianMap = 1 - exp(-((X - sizeX/3).^2 + (Y - sizeY/2).^2)/(2*10^2));
interpolateGaussian1 = @(x,y) interp2((2:1:sizeX-1),(2:1:sizeY-1),gaussianMap,x,y);

gaussianMap2 = 1 - exp(-((X - 30).^2 + (Y - 50).^2)/(2*15^2));
interpolateGaussian2 = @(x,y) interp2((2:1:sizeX-1),(2:1:sizeY-1),gaussianMap2,x,y);

iif = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();
flowX = @(x,y) (cos(atan2(y - positionY(1),x - 45) + pi/2) ./ bsxfun(@max, ones(size(x)),abs(((x - 45).^2 + (y - positionY(1)).^2).^0.5 - radius(1))).^0.45) ...
               .* interpolateGaussian2(x, y);%.* interpolateGaussian2(x, y);
flowY = @(x,y) (sin(atan2(y - positionY(1),x - 45) + pi/2) ./ bsxfun(@max, ones(size(x)),abs(((x - 45).^2 + (y - positionY(1)).^2).^0.5 - radius(1))).^0.45) ...
               .* interpolateGaussian2(x, y);%.* interpolateGaussian2(x, y);

noiseSize = 1;
driftSize = 5;
gradientSize = 15;

figureHandle = figure(1);
figureHandle.Renderer='Painters';
clf;
imagesc(energyLandscape);
colormap(gray)
hold on;

x = (2:3:sizeX-1);
y = (2:3:sizeY-1);
[X, Y] = meshgrid(x, y);
X = reshape(X, 1, size(X,1)*size(X,2));
Y = reshape(Y, 1, size(X,1)*size(X,2));
arrowX = reshape(gradUx(X, Y), 1, size(X,1)*size(X,2));
arrowY = reshape(gradUy(X, Y), 1, size(X,1)*size(X,2));
% arrowFlowX = reshape(flowX(sub2ind(size(flowX), X, Y)), 1, size(X,1)*size(X,2));
% arrowFlowY = reshape(flowY(sub2ind(size(flowY), X, Y)), 1, size(X,1)*size(X,2));
arrowFlowX = reshape(flowX(X, Y), 1, size(X,1)*size(X,2));
arrowFlowY = reshape(flowY(X, Y), 1, size(X,1)*size(X,2));
% quiver(X, Y, arrowX, arrowY, 1);
% quiver(X, Y, arrowFlowX, arrowFlowY, 3, 'r');
colormap jet;

validTrace = false;
while ~validTrace
    trace = [];
%     x = floor(rand(1,2) .* [sizeX, sizeY]);
    x = [20, 50];
    for i = 1:1000
        lastX = x;

        driftX = -gradientSize * gradUx(x(1), x(2)) + driftSize * flowX(x(1), x(2));
        driftY = -gradientSize * gradUy(x(1), x(2)) + driftSize * flowY(x(1), x(2));
        noise = noiseSize * normrnd(0, 1, 1, 2);

        x = x + [driftX, driftY] + noise;

        %x = x - [gradUx(roundedX(1), roundedX(2)), gradUy(roundedX(1), roundedX(2))] * potentialSize + normrnd(0, 2, 1, 2) * noiseSize;
        %x = x - [flowX(roundedX(1), roundedX(2)), flowY(roundedX(1), roundedX(2))] * potentialSize + normrnd(0, 2, 1, 2) * noiseSize;

        x(1) = max(min(x(1), sizeX-1), 2);
        x(2) = max(min(x(2), sizeY-1), 2);

        line = [lastX; x];
        plot(line(:,1), line(:,2), 'k-');
    %     
    %     pause(0.0005);

        trace(:,i) = x;
    end

    validTrace = true;
%     leftSides = find(trace(1,300:end) < positionX(1));
%     rightSides = find(trace(1,300:end) > positionX(3));
%     if range(leftSides) < 50 || range(rightSides) < 50
%         validTrace = false;
%     end
end

% t = 1:300;
% period = 120;
% startRadius = 10;
% incRadius = 10;
% radiusSigma = 0;
% rollX = @(t) (startRadius + incRadius*t/period).*cos(t*2*pi/period) + normrnd(0, radiusSigma, size(t));
% rollY = @(t) (startRadius + incRadius*t/period).*sin(t*2*pi/period) + normrnd(0, radiusSigma, size(t));
% plot(rollX(t), rollY(t));
% 
finalTrace = trace;
% finalTrace = [rollX(trace(1,:)); rollY(trace(1,:)); trace(2,:)];

figureHandle = figure(2);
figureHandle.Renderer='Painters';
clf
hold on
plot(finalTrace(1,:));
plot(finalTrace(2,:));

%% Different trough
seed = 7;

startPoints = [49 30; 49 24];
trialNumber = 10;

trace = [];
for j = 1:size(startPoints,1)
    for k = 1:trialNumber
        rng(seed * j * k);
        
        x = startPoints(j,:);
        for i = 1:60
            lastX = x;

            driftX = -gradientSize * gradUx(x(1), x(2)) + driftSize * flowX(x(1), x(2));
            driftY = -gradientSize * gradUy(x(1), x(2)) + driftSize * flowY(x(1), x(2));
            noise = noiseSize * normrnd(0, 1, 1, 2);

            x = x + [driftX, driftY] + noise;

            x(1) = max(min(x(1), sizeX-1), 2);
            x(2) = max(min(x(2), sizeY-1), 2);

            trace(:,i,j,k) = x;
        end
    end
end

figureHandle = figure(3);
figureHandle.Renderer='Painters';
clf
contour(-energyLandscape);
hold on
for i = 1:size(trace,4)
    plot(trace(1,:,1,i), trace(2,:,1,i), 'Color', 'r');
end
for i = 1:size(trace,4)
    plot(trace(1,:,2,i), trace(2,:,2,i), 'Color', 'b');
end

%% Same trough
seed = 7;

startPoints = [64 30; 70 24];
trialNumber = 10;

trace = [];
for j = 1:size(startPoints,1)
    for k = 1:trialNumber
        rng(seed * j * k);
        
        x = startPoints(j,:);
        for i = 1:50
            lastX = x;

            driftX = -gradientSize * gradUx(x(1), x(2)) + driftSize * flowX(x(1), x(2));
            driftY = -gradientSize * gradUy(x(1), x(2)) + driftSize * flowY(x(1), x(2));
            noise = noiseSize * normrnd(0, 1, 1, 2);

            x = x + [driftX, driftY] + noise;

            x(1) = max(min(x(1), sizeX-1), 2);
            x(2) = max(min(x(2), sizeY-1), 2);

            trace(:,i,j,k) = x;
        end
    end
end

figureHandle = figure(4);
figureHandle.Renderer='Painters';
clf
contour(-energyLandscape);
hold on
for i = 1:size(trace,4)
    plot(trace(1,:,1,i), trace(2,:,1,i), 'Color', 'r');
end
for i = 1:size(trace,4)
    plot(trace(1,:,2,i), trace(2,:,2,i), 'Color', 'b');
end

%% Same trough no noise
seed = 7;

startPoints = [64 30; 70 24];
trialNumber = 10;

trace = [];
for j = 1:size(startPoints,1)
    for k = 1:trialNumber
        rng(seed * j * k);
        
        x = startPoints(j,:);
        for i = 1:30
            lastX = x;

            driftX = -gradientSize * gradUx(x(1), x(2)) + driftSize * flowX(x(1), x(2));
            driftY = -gradientSize * gradUy(x(1), x(2)) + driftSize * flowY(x(1), x(2));
            noise = noiseSize * normrnd(0, 1, 1, 2);

            if i == 1
                x = x + [driftX, driftY] + noise;
            else
                x = x + [driftX, driftY];
            end

            x(1) = max(min(x(1), sizeX-1), 2);
            x(2) = max(min(x(2), sizeY-1), 2);

            trace(:,i,j,k) = x;
        end
    end
end

figureHandle = figure(5);
figureHandle.Renderer='Painters';
clf
contour(-energyLandscape);
hold on
for i = 1:size(trace,4)
    plot(trace(1,:,1,i), trace(2,:,1,i), 'Color', 'r');
end
for i = 1:size(trace,4)
    plot(trace(1,:,2,i), trace(2,:,2,i), 'Color', 'b');
end


%%

noise = 0.07;
restoreForce = 0.03;

starts = [-2 -1 0 1 2];

colors = parula(length(starts));


traces = starts;
for t = 2:100
    traces(t,:) = traces(t-1,:) - restoreForce*sign(traces(t-1,:)).*(abs(traces(t-1,:)))  + noise*normrnd(0,1,size(starts));
end


figureHandle = figure(3);
figureHandle.Renderer='Painters';
clf
hold on
for i = 1:length(starts)
    plot(traces(:,i), 'Color', colors(i,:));
end


traces = starts;
for t = 2:100
    traces(t,:) = traces(t-1,:) - restoreForce*sign(traces(t-1,:)).*(abs(traces(t-1,:)))  + 0.00*normrnd(0,1,size(starts));
end


figureHandle = figure(4);
figureHandle.Renderer='Painters';
clf
hold on
for i = 1:length(starts)
    plot(traces(:,i), 'Color', colors(i,:));
end


targets = starts/2;
traces = starts;
for t = 2:100
    traces(t,:) = traces(t-1,:) - restoreForce*sign(traces(t-1,:) - targets).*(abs(traces(t-1,:) - targets))  + noise*normrnd(0,1,size(starts));
end


figureHandle = figure(5);
figureHandle.Renderer='Painters';
clf
hold on
plot(traces);
for i = 1:length(starts)
    plot(traces(:,i), 'Color', colors(i,:));
end


%% Display results

PLOT_DELAY_EMBEDDED = 0;

app.SavedData = saveData;

numTrial = max(app.SavedData.TrialData);
trialLength = size(app.SavedData.FinalStream,1) / numTrial;
times = 1:trialLength;

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

figureHandle = figure(7);
figureHandle.Renderer='Painters';
clf;
hold on;

colors = jet(7);
colors = colors([1 2 4 5 6 7],:);
colors(3,:) = colors(3,:)*0.6;
colors(4,:) = [238 210 2]/256;
% colors(4,:) = colors(3,:);

lineColors = [4 2];
lineColors = colors(lineColors,:);

% colors = jet(app.SavedData.BestLoopCount);
% h = plot3(finalStream(loopStarts+1,:)*trialPCABasis(:,1), finalStream(loopStarts+1,:)*trialPCABasis(:,2), finalStream(loopStarts+1,:)*trialPCABasis(:,3), 'o', 'LineWidth', 0.5);

% for i = [6, 7]%app.SavedData.BestLoopCount    
% for i = [4, 5]%app.SavedData.BestLoopCount    
% for i = [1,2,5,6]%app.SavedData.BestLoopCount    
% for i = [3,7,8,10]%app.SavedData.BestLoopCount    
allClusterMeans = [];
for i = 1:app.SavedData.BestLoopCount
    thisLoopIDs = find(app.SavedData.BestLoopAssignments(:,1) == i);
    thisLoopClusters = app.SavedData.BestLoopAssignments(thisLoopIDs,2);
    
    allIDs = find(app.SavedData.BestStateMap(:,1) == i);
    trialIDs = floor((allIDs-1)/trialLength) + 1;
    conditionID = i;
    
    badIDs = find(app.SavedData.BestStateMap(:,1) == i);
    badIDs = setdiff(1:size(finalStream,1), badIDs);
    
    plotStream = finalStream;
    plotStream(badIDs,:) = nan;
    
    if PLOT_DELAY_EMBEDDED
        h = plot3(plotStream(trialIndicies,:)*trialPCABasis(:,1), plotStream(trialIndicies,:)*trialPCABasis(:,2), plotStream(trialIndicies,:)*trialPCABasis(:,3), 'LineWidth', 0.5, 'Color', lineColors(conditionID,:));
    else
        h = plot(plotStream(trialIndicies,1), plotStream(trialIndicies,2), 'LineWidth', 0.5, 'Color', lineColors(conditionID,:));
    end
    h.Color(4) = 0.9;
    
    meanTimes = [];
    clusterLengths = [];
    clusterMeans = [];
    clusterSTDs = [];
    for j = 1:length(thisLoopClusters)
        thisIndices = find(app.SavedData.BestStateMap(:,1) == i & app.SavedData.BestStateMap(:,2) == thisLoopClusters(j));

        meanTimes(j) = mode(mod(thisIndices, trialLength)+1);
        clusterLengths(j) = length(thisIndices);
        
        clusterMeans(j,:) = nanmean(plotStream(thisIndices, :), 1);
        clusterSTDs(j,:) = nanstd(plotStream(thisIndices, :), [], 1);
        
        allClusterMeans(i, thisLoopClusters(j),:) = clusterMeans(j,:);
    end
    
%     startPoints = find(app.SavedData.BestStateMap(loopStarts,1) == i);
%     bestStartPoint = mode(app.SavedData.BestStateMap(loopStarts(startPoints),2));
%     
%     
    thisLoopIDs = 1:length(thisLoopClusters);
    clusterOrder = 1:length(thisLoopClusters);
%     clusterOrder = app.SavedData.BestLoopAssignments(thisLoopIDs,2);
    
%     testTimes = meanTimes;
%     testTimes(clusterLengths < 5) = 1000;
%     [~, startCluster] = min(meanTimes);
    startCluster = 1;
    
    sortedOrder = [startCluster:length(clusterOrder) 1:startCluster-1];
    thisLoopIDs = thisLoopIDs(sortedOrder);
    
%     for j = 1:length(clusterOrder)
%         thisIndices = find(app.SavedData.BestLoopAssignments(:,1) == i & app.SavedData.BestLoopAssignments(:,2) == clusterOrder(sortedOrder(j)));
%     end
    
    meanTimes = meanTimes(sortedOrder);
    
    badTimes = find(meanTimes <= min(times) | meanTimes >= max(times));
    goodTimes = setdiff(sortedOrder, badTimes);
    
%     thisTrace = app.SavedData.BestEmission(thisLoopIDs(goodTimes), :);
    clusterSTDs = clusterSTDs(thisLoopIDs(goodTimes),:);
    thisTrace = clusterMeans(thisLoopIDs(goodTimes),:);
    thisTrace = filterData(thisTrace', 0.5, [], 1, 0)';
    
    if PLOT_DELAY_EMBEDDED
%         plot3(thisTrace*trialPCABasis(:,1), thisTrace*trialPCABasis(:,2), thisTrace*trialPCABasis(:,3), 'LineWidth', 2, 'Color', colors(i,:))
    
        tubeMean = thisTrace*trialPCABasis(:,1:3);
        projectedSTDs = clusterSTDs*trialPCABasis(:,1:3);
    else
%         plot(thisTrace(:,1), thisTrace(:,2), 'LineWidth', 2, 'Color', colors(i,:))
        
        tubeMean = [thisTrace(:,1:2) zeros(size(thisTrace,1), 1)];
        projectedSTDs = [clusterSTDs(:,1:2) zeros(size(thisTrace,1), 1)];
    end
    
    projectedDerivative = diff(tubeMean);
    projectedDerivative = [projectedDerivative(1,:); projectedDerivative];
    
    orthogonalSTDs = zeros(size(projectedSTDs,1),1);
    for j = 1:size(projectedDerivative,1)
%         normDerivative = projectedDerivative(j,:)/norm(projectedDerivative(j,:));
%         derivativeNullspace = null(normDerivative)';
%         
%         changeOfBasis = [normDerivative; derivativeNullspace];
%         changedSTD = changeOfBasis * projectedSTDs(j,:)';
%         
%         orthogonalSTDs(j) = norm(changedSTD(2:end));

        orthogonalSTDs(j) = norm(projectedSTDs(j,:));
    end
    
    orthogonalSTDs = filterData(orthogonalSTDs, 0.5, [], 1, 0);
    
    [X,Y,Z,V] = tubeplot(tubeMean(:,1), tubeMean(:,2), tubeMean(:,3), orthogonalSTDs, ones(1,size(tubeMean,1)),10,[0 0 1]);

    t = surf(X,Y,Z,V, 'FaceColor', lineColors(conditionID,:));
    if conditionID > 3
        t.EdgeColor = [0 0 0];
        t.EdgeAlpha = 0.7;
    else
        t.EdgeColor = lineColors(conditionID,:);
        t.EdgeAlpha = 0.5;
    end
    t.FaceAlpha = 0.2;
end

%% Display results phase

PLOT_DELAY_EMBEDDED = 0;

app.SavedData = saveData;

numTrial = max(app.SavedData.TrialData);
trialLength = size(app.SavedData.FinalStream,1) / numTrial;
times = 1:trialLength;

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

figureHandle = figure(7);
figureHandle.Renderer='Painters';
clf;
hold on;

colors = jet(7);
colors = colors([1 2 4 5 6 7],:);
colors(3,:) = colors(3,:)*0.6;
colors(4,:) = [238 210 2]/256;
% colors(4,:) = colors(3,:);

lineColors = [4 2];
lineColors = colors(lineColors,:);

% colors = jet(app.SavedData.BestLoopCount);
% h = plot3(finalStream(loopStarts+1,:)*trialPCABasis(:,1), finalStream(loopStarts+1,:)*trialPCABasis(:,2), finalStream(loopStarts+1,:)*trialPCABasis(:,3), 'o', 'LineWidth', 0.5);

% for i = [6, 7]%app.SavedData.BestLoopCount    
% for i = [4, 5]%app.SavedData.BestLoopCount    
% for i = [1,2,5,6]%app.SavedData.BestLoopCount    
% for i = [3,7,8,10]%app.SavedData.BestLoopCount    
allClusterMeans = [];
for i = 1:app.SavedData.BestLoopCount
    thisLoopIDs = find(app.SavedData.BestLoopAssignments(:,1) == i);
    thisLoopClusters = app.SavedData.BestLoopAssignments(thisLoopIDs,2);
    
    allIDs = find(app.SavedData.BestStateMap(:,1) == i);
    trialIDs = floor((allIDs-1)/trialLength) + 1;
    conditionID = i;
    
    badIDs = find(app.SavedData.BestStateMap(:,1) == i);
    badIDs = setdiff(1:size(finalStream,1), badIDs);
    
    plotStream = finalStream;
    plotStream(badIDs,:) = nan;
    
%     if PLOT_DELAY_EMBEDDED
%         h = plot3(plotStream(trialIndicies,:)*trialPCABasis(:,1), plotStream(trialIndicies,:)*trialPCABasis(:,2), plotStream(trialIndicies,:)*trialPCABasis(:,3), 'LineWidth', 0.5, 'Color', lineColors(conditionID,:));
%     else
%         h = plot(plotStream(trialIndicies,1), plotStream(trialIndicies,2), 'LineWidth', 0.5, 'Color', lineColors(conditionID,:));
%     end
%     h.Color(4) = 0.9;
    
    meanTimes = [];
    clusterLengths = [];
    clusterMeans = [];
    clusterSTDs = [];
    for j = 1:length(thisLoopClusters)
        thisIndices = find(app.SavedData.BestStateMap(:,1) == i & app.SavedData.BestStateMap(:,2) == thisLoopClusters(j));

        meanTimes(j) = mode(mod(thisIndices, trialLength)+1);
        clusterLengths(j) = length(thisIndices);
        
        clusterMeans(j,:) = nanmean(plotStream(thisIndices, :), 1);
        clusterSTDs(j,:) = nanstd(plotStream(thisIndices, :), [], 1);
        
        allClusterMeans(i, thisLoopClusters(j),:) = clusterMeans(j,:);
    end
    
%     startPoints = find(app.SavedData.BestStateMap(loopStarts,1) == i);
%     bestStartPoint = mode(app.SavedData.BestStateMap(loopStarts(startPoints),2));
%     
%     
    thisLoopIDs = 1:length(thisLoopClusters);
    clusterOrder = 1:length(thisLoopClusters);
%     clusterOrder = app.SavedData.BestLoopAssignments(thisLoopIDs,2);
    
%     testTimes = meanTimes;
%     testTimes(clusterLengths < 5) = 1000;
%     [~, startCluster] = min(meanTimes);
    startCluster = 1;
    
    sortedOrder = [startCluster:length(clusterOrder) 1:startCluster-1];
    thisLoopIDs = thisLoopIDs(sortedOrder);
    
%     for j = 1:length(clusterOrder)
%         thisIndices = find(app.SavedData.BestLoopAssignments(:,1) == i & app.SavedData.BestLoopAssignments(:,2) == clusterOrder(sortedOrder(j)));
%     end
    
    meanTimes = meanTimes(sortedOrder);
    
    badTimes = find(meanTimes <= min(times) | meanTimes >= max(times));
    goodTimes = setdiff(sortedOrder, badTimes);
    
%     thisTrace = app.SavedData.BestEmission(thisLoopIDs(goodTimes), :);
    clusterSTDs = clusterSTDs(thisLoopIDs(goodTimes),:);
    thisTrace = clusterMeans(thisLoopIDs(goodTimes),:);
    thisTrace = filterData(thisTrace', 0.5, [], 1, 0)';
    
    if PLOT_DELAY_EMBEDDED
%         plot3(thisTrace*trialPCABasis(:,1), thisTrace*trialPCABasis(:,2), thisTrace*trialPCABasis(:,3), 'LineWidth', 2, 'Color', colors(i,:))
    
        tubeMean = thisTrace*trialPCABasis(:,1:3);
        projectedSTDs = clusterSTDs*trialPCABasis(:,1:3);
    else
%         plot(thisTrace(:,1), thisTrace(:,2), 'LineWidth', 2, 'Color', colors(i,:))
        
        tubeMean = [thisTrace(:,1:2) zeros(size(thisTrace,1), 1)];
        projectedSTDs = [clusterSTDs(:,1:2) zeros(size(thisTrace,1), 1)];
    end
    
    projectedDerivative = diff(tubeMean);
    projectedDerivative = [projectedDerivative(1,:); projectedDerivative];
    
    orthogonalSTDs = zeros(size(projectedSTDs,1),1);
    for j = 1:size(projectedDerivative,1)
%         normDerivative = projectedDerivative(j,:)/norm(projectedDerivative(j,:));
%         derivativeNullspace = null(normDerivative)';
%         
%         changeOfBasis = [normDerivative; derivativeNullspace];
%         changedSTD = changeOfBasis * projectedSTDs(j,:)';
%         
%         orthogonalSTDs(j) = norm(changedSTD(2:end));

        orthogonalSTDs(j) = norm(projectedSTDs(j,:));
    end
    
    orthogonalSTDs = filterData(orthogonalSTDs, 0.5, [], 1, 0);
    
    midPosition = [40 50];
    phases = atan2(tubeMean(:,2) - midPosition(2), tubeMean(:,1) - midPosition(1));
    
    [X,Y,Z,V] = tubeplot(tubeMean(:,1), tubeMean(:,2), tubeMean(:,3), orthogonalSTDs, ones(1,size(tubeMean,1)),10,[0 0 1]);

    colors = hsv(256);
    discreteTimes = floor((phases + pi) / (2*pi) * 256) + 1;
    colors = reshape(colors(discreteTimes,:), [size(X,1) 1 3]);
    colors = repmat(colors, [1 size(X,2), 1]);

    t = surf(X,Y,Z,colors);
    t.EdgeAlpha = 0.2;
    t.FaceAlpha = 0.5;
    
    xlim([0 90])
    ylim([10 100])
end

%% Markov matrix

figureHandle = figure(8);
figureHandle.Renderer='Painters';
clf;
colormap('gray');
imagesc(1 - saveData.BestModel);
hold on;

phaseValues = [];
for i = 1:size(saveData.BestLoopAssignments)
    if saveData.BestLoopAssignments(i,1) == 1
        center = [40 50 0];
    else
        center = [55 55 0];
    end
    
    thisMean = squeeze(allClusterMeans(saveData.BestLoopAssignments(i,1), saveData.BestLoopAssignments(i,2),:))'*trialPCABasis(:,1:3);;
    
    thisDiff = thisMean - center;
    
    thisAngle = atan2(thisDiff(2), thisDiff(1));
    
    phaseValues(i) = thisAngle;
end

figureHandle = figure(9);
figureHandle.Renderer='Painters';
clf;
subplot(1,2,1)
plot(saveData.BestLoopAssignments(:,1));
subplot(1,2,2)
plot(phaseValues);
hold on;

figureHandle = figure(10);
figureHandle.Renderer='Painters';
clf;
colormap(hsv)
imagesc(phaseValues);
hold on;



