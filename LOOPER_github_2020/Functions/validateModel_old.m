
% finalDynamicsStream
% validationModel
% validationEmission

% returns likelihood (bootstrapped likelihood

if size(validationEmission, 2) ~= size(finalDynamicsStream,2)
    disp(['Emission matrix must be STATES by CHANNELS (' num2str(size(validationModel,2)) ', ' num2str(size(finalDynamicsStream,2)) ')']);
    return;
end

if size(validationEmission, 1) ~= size(validationModel,2)
    disp(['Emission matrix must be STATES by CHANNELS (' num2str(size(validationModel,2)) ', ' num2str(size(finalDynamicsStream,2)) ')']);
    return;
end

disp('Validating...');

BOOTSTRAP_AMOUNT = 5;
USE_PCA = 1;

likelihoods = zeros(1, BOOTSTRAP_AMOUNT);

noInputDynamics = finalDynamicsStream;%(:,1:100);
noInputEmissions = validationEmission;%(:,1:100);

[finalPCABasis, ~, ~, ~, explained] = pca(noInputDynamics, 'centered', false);
requiredDimensions = round(find(cumsum(explained) > 95, 1, 'first'));
if USE_PCA
    finalPCABasis = finalPCABasis(:, 1:requiredDimensions);
else
    finalPCABasis = eye(size(finalPCABasis, 1));
end

for bootstrapIndex = 1:BOOTSTRAP_AMOUNT
    SIMULATION_LENGTH = size(noInputDynamics,1)*10;
    simulatedStates = [];
    simulatedStates(1) = 1;
    for i = 2:SIMULATION_LENGTH
        simulatedStates(i) = randsample(1:size(validationModel,1), 1, true, validationModel(simulatedStates(i-1),:));
    end

    simulatedData = noInputEmissions(simulatedStates,:);
%     simulatedData = filterData(simulatedData, 1);

    normalizedFrequencies = 2*pi ./ logspace(1, log10(size(noInputDynamics,1)/2), 30);
    observedSpectrum = [];
    simulatedSpectrum = [];
    observedSpectrumCI = [];
    simulatedSpectrumCI = [];
    for i = 1:size(finalPCABasis, 2)
        [observedSpectrum(i,:),~,observedSpectrumCI(i,:,:)] = pmtm(noInputDynamics*finalPCABasis(:,i), 10, normalizedFrequencies,'ConfidenceLevel',0.68);
        [simulatedSpectrum(i,:),~,simulatedSpectrumCI(i,:,:)] = pmtm(simulatedData*finalPCABasis(:,i), 10, normalizedFrequencies,'ConfidenceLevel',0.68);
    end

    observedSTD = observedSpectrum - observedSpectrumCI(:,:,2);
    simulatedSTD = simulatedSpectrum - simulatedSpectrumCI(:,:,2);

    pValues = [];
    for i = 1:size(finalPCABasis, 2)
        differences = abs(observedSpectrum(i,:) - simulatedSpectrum(i,:));
        zscores = differences ./ sqrt(observedSTD(i,:).^2 + simulatedSTD(i,:).^2);
        pValues(i,:) = (1 - normcdf(zscores));
    end

    likelihood = -mean(log(pValues),2);
%     likelihood
    
    if USE_PCA
        likelihoodWeights = explained(1:length(likelihood));
    else
        likelihoodWeights = ones(length(likelihood),1);
    end
    likelihoodWeights = likelihoodWeights / sum(likelihoodWeights);
    likelihood = sum(likelihood);% .* likelihoodWeights);

    likelihoods(bootstrapIndex) = likelihood;
end

figure(1);
clf;
clear h;
h(2) = subplot(2,1,2);
imagesc((simulatedData(1:size(noInputDynamics,1),:)*finalPCABasis)');
title('Simulated');
colorAxis = caxis;
h(1) = subplot(2,1,1);
imagesc((noInputDynamics*finalPCABasis)');
caxis(colorAxis);
title('Observed');
linkaxes(h, 'xy');

plotChannel = 1;

for i = 1:10
    plotChannel = i;
    
    figure(i+1);
    clf;
    hold on
    plot(normalizedFrequencies, observedSpectrum(plotChannel,:), 'b', 'lineWidth', 2);
    plot(normalizedFrequencies, observedSpectrumCI(plotChannel,:,1), 'b');
    plot(normalizedFrequencies, observedSpectrumCI(plotChannel,:,2), 'b');
    plot(normalizedFrequencies, simulatedSpectrum(plotChannel,:), 'r', 'lineWidth', 2);
    plot(normalizedFrequencies, simulatedSpectrumCI(plotChannel,:,1), 'r');
    plot(normalizedFrequencies, simulatedSpectrumCI(plotChannel,:,2), 'r');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
end