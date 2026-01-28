
dataDirectory = 'F:/Dropbox/FullBrain/wbdata/';

load([dataDirectory 'neuronList18']);

allFiles = dir([dataDirectory 'TS*.mat']);

allData = {};
for i = 1:length(allFiles)
    load([dataDirectory allFiles(i).name]);
    
    rawData = wbData.deltaFOverF_bc';
    
    rawData = detrend(rawData);
    
    endIndex = strfind(wbData.FlNm, '_');
    endIndex = endIndex(1);
    
    switch (wbData.FlNm(endIndex-1))
        case 'b'
            wormID = 5;

        case 'c'
            wormID = 3;

        case 'd'
            wormID = 4;

        case 'e'
            wormID = 1;

        case 'f'
            wormID = 2;
    end
    
    sharedData = rawData(neuronList(wormID(1),:),:);
    sharedData = zscore(sharedData, 0 ,2);
    
    allData{i} = sharedData;
end

concatData = [];
for i = 1:length(allData)
    concatData = [concatData allData{i}];
end
[pcaBasis, ~, ~, ~, explained] = pca(concatData);