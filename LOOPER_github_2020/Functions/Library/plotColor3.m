
function plotColor3(xValues, yValues, zValues, colorIDs, colorMap, lineWidth, alpha)

IDs = unique(colorIDs);

if ~exist('colorMap', 'var') || isempty(colorMap)
    colorMap = parula(length(IDs));
end

if ~exist('lineWidth', 'var') || isempty(lineWidth)
    lineWidth = 3;
end

if ~exist('alpha', 'var') || isempty(alpha)
    alpha = 1;
end

hold on;
for j = 1:length(IDs)
    ID = IDs(j);
    indices = find(colorIDs == ID);
    startIndex = 0;
    endIndex = 0;
    for i = 1:length(indices)-1
        if startIndex == 0
            startIndex = indices(i);
        end
        if indices(i + 1) > indices(i) + 1 || i == length(indices)-1
            endIndex = indices(i);
        end

        if endIndex ~= 0
            trace = [startIndex-1:endIndex+1];
            trace(trace <= 0) = [];
            p = plot3(xValues(trace), yValues(trace), zValues(trace), 'LineWidth',lineWidth, 'Color', colorMap(mod(ID-1,size(colorMap,1))+1,:));
            p.Color(4) = alpha;
            
            startIndex = 0;
            endIndex = 0;
        end
    end
end
