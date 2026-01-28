
% returns scoreMean, scoreSTD

maxTime = 10;

scoreMean = 0;
scoreSTD = 0;
FIG_BASE = 300;

if size(validationEmission, 2) ~= size(finalDynamicsStream,2)
    disp(['Test data must have ' num2str(originalDataSize(1)) ' CHANNELS']);
    return;
end

if size(validationEmission, 1) ~= size(validationModel,2)
    disp(['Emission matrix must be STATES by CHANNELS (' num2str(size(validationModel,2)) ', ' num2str(size(finalDynamicsStream,1)) ')']);
    return;
end

disp('Validating...');

USE_PCA = 0;

noInputDynamics = finalDynamicsStream;
noInputEmissions = validationEmission;

dynamicsMean = mean(noInputDynamics,1);
dynamicsSTD = std(noInputDynamics - dynamicsMean,1);

noInputDynamics = (noInputDynamics - dynamicsMean) ./ dynamicsSTD;
noInputEmissions = (noInputEmissions - dynamicsMean) ./ dynamicsSTD;

noInputDynamics(isnan(noInputDynamics)) = 0;
noInputEmissions(isnan(noInputEmissions)) = 0;

if USE_PCA
    [finalPCABasis, ~, ~, ~, explained] = pca(noInputDynamics);
    requiredDimensions = round(find(cumsum(explained) > 95, 1, 'first'));

    finalPCABasis = finalPCABasis(:, 1:requiredDimensions);
else
    finalPCABasis = eye(size(noInputDynamics, 2));
end

emissionDistances = pdist2(noInputEmissions, noInputDynamics).^2;
emissionVelocities = repmat(permute(noInputEmissions, [3,1,2]), [size(noInputEmissions,1),1,1]) - repmat(permute(noInputEmissions, [1,3,2]), [1,size(noInputEmissions,1),1]);

predictTime = 1;

possibleTransitions = validationModel;
possibleTransitions = possibleTransitions^predictTime;

fromStates = [];
toStates = [];
scores = [];
for i = 1:size(noInputDynamics,1)-predictTime
    distances = emissionDistances(:,i);
    nextDistances = emissionDistances(:,i+predictTime);
    
    transitionScores = log(((distances .* repmat(nextDistances, [1, size(nextDistances,1)])) ./ possibleTransitions));
    transitionScores(isinf(transitionScores)) = 100000;
    
    scores(i) = nanmin(transitionScores(:));
    
    index = find(transitionScores == scores(i), 1);
    [currentState, nextState] = ind2sub(size(transitionScores), index);
    
    fromStates(i) = currentState;
    toStates(i) = nextState;
end


scoreMean = mean(scores);
scoreSTD = std(scores);

%%

if size(finalDynamicsStream,2) > 2
    [pcaBasis, pcaOutputs] = pca(finalDynamicsStream, 'NumComponents', 3);
else
    pcaBasis = eye(size(finalDynamicsStream,2));
end

colors = jet(256);
scoreIndices = scores - min(scores);
scoreIndices = scoreIndices / max(scoreIndices);
scoreIndices = round(1 + scoreIndices * 255);

figure(FIG_BASE + 1);
clf;
hold on;
if size(pcaBasis, 2) > 2
    h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
else
    h = plot(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), 'LineWidth', 0.5);
end
h.Color(4) = 0.2;

plotDynamics = finalDynamicsStream(1:end-predictTime,:);

if size(pcaBasis, 2) > 2
    scatter3(plotDynamics*pcaBasis(:,1), plotDynamics*pcaBasis(:,2), plotDynamics*pcaBasis(:,3), 32, scores);
else
    scatter(plotDynamics*pcaBasis(:,1), plotDynamics*pcaBasis(:,2), 32, scores);
end
colormap(jet)

if size(pcaBasis, 2) > 2
    scatter3(validationEmission*pcaBasis(:,1), validationEmission*pcaBasis(:,2), validationEmission*pcaBasis(:,3), 100, 'kx', 'LineWidth', 3);
else
    scatter(validationEmission*pcaBasis(:,1), validationEmission*pcaBasis(:,2), 100, 'kx', 'LineWidth', 3);
end

for i = 1:5:length(toStates)
    maxDims = min(3, size(pcaBasis, 2));
    lines = validationEmission([fromStates(i),toStates(i)],:)*pcaBasis(:,1:maxDims);
    
    if size(pcaBasis, 2) > 2
        plot3(lines(:,1), lines(:,2), lines(:,3), 'k', 'LineWidth', 1);
    else
        plot(lines(:,1), lines(:,2), 'k', 'LineWidth', 1);
    end
end

title('Validation: dynamics, emissions, and best transitions');
if size(pcaBasis, 2) <= 2
    xlabel('PC1');
    ylabel('PC2');
end
