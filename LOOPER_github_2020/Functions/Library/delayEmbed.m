function [delayEmbeddedData, embeddedTime] = delayEmbed(data, numDelays, deltaT, useDerivatives, dontFlip)
    if ~exist('useDerivatives') || isempty(useDerivatives)
        useDerivatives = 1;
    end
    
    if ~exist('dontFlip') || isempty(dontFlip)
        dontFlip = 0;
    end
    
    if numDelays == 0 || deltaT == 0
        delayEmbeddedData = data;
        embeddedTime = 0;
        return;
    end
    
    flipped = 0;
    if ~dontFlip
        if size(data,2) < size(data,1)
            data = data';
            flipped = 1;
        end
    end

    dataDiffs = [];
    if useDerivatives
        filteredTrace = diff(data,1,2);
        
        dataDiffs(:,:) = [zeros(size(data,1),1), filteredTrace];
    end

    %%
    embeddedTime = numDelays*deltaT + useDerivatives;
    index = 1;
    for i = 1 + embeddedTime:size(data, 2)
        thisVector = [];

        for j = 0:numDelays
            thisVector = [thisVector; data(:,i - j*deltaT)];

            if useDerivatives
                thisVector = [thisVector; dataDiffs(:,i - j*deltaT)];
            end
        end

        delayEmbeddedData(:,index) = thisVector;
        index = index + 1;
    end
    
    if flipped
        delayEmbeddedData = delayEmbeddedData';
    end