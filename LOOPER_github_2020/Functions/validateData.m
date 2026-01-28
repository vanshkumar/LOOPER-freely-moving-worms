function [] = validateData(LOOPERData, validateData, validateInputs, validateOutputs)
    rawData = validateData;
    rawData = convertToCell(rawData);

    lastEnd = 0;
    for i = 1:length(rawData)
        trialData(lastEnd + (1:size(rawData{i}, 2))) = i;

        lastEnd = lastEnd + size(rawData{i}, 2);
    end

    rawData = mergeData(rawData);
    
    inputs = [];
    if exist('validateInputs') && ~isempty(validateInputs)
        inputs = mergeData(convertToCell(validateInputs));
    end
    
    outputs = [];
    if exist('validateOutputs') && ~isempty(validateOutputs)
        outputs = mergeData(convertToCell(validateOutputs));
    end
    
    [altDynamics, ~, ~, ~, ~] = preprocessData(rawData, inputs, LOOPERData.PreprocessData.InputLambda, outputs, LOOPERData.PreprocessData.OutputLambda, trialData, LOOPERData.PreprocessData.TrialLambda, false, LOOPERData.PreprocessData.Smoothing, LOOPERData.PreprocessData.ZScore, LOOPERData.PreprocessData.DelayTime, LOOPERData.PreprocessData.DelayCount, LOOPERData.DataMean, LOOPERData.DataSTD);

    finalDynamicsStream = altDynamics';
    validationModel = LOOPERData.BestModel;
    validationEmission = LOOPERData.BestEmission;

    originalDataSize = size(LOOPERData.RawData);

    validateModel;
    
    numDataPoints = size(finalDynamicsStream, 1);
            
    disp(['Score: ' num2str(scoreMean) ' +/- ' num2str(1.96 * scoreSTD / sqrt(size(finalDynamicsStream,1)))]);
end
