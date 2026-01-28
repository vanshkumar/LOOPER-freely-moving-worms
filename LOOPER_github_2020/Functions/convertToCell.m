function [cellData] = convertToCell(data)
    if size(data,3) > 1
        cellData = cell(size(data,3),1);
        for i = 1:size(data,3)
            cellData{i} = data(:,:,i);
        end
    else
        if ~iscell(data)
            cellData{1} = data;
        else
            cellData = data;
        end
    end
end