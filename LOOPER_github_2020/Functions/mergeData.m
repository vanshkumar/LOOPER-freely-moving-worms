function [mergedData] = mergeData(data)
    totalLength = 0;

    for i = 1:length(data)
        totalLength = totalLength + size(data{i},2);
    end

    mergedData = zeros(size(data{1}, 1), totalLength);
    currentEnd = 0;
    for i = 1:length(data)
        mergedData(:,(1:size(data{i},2))+currentEnd) = data{i};
        currentEnd = currentEnd + size(data{i},2);
    end
end