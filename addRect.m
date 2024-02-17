% places a rectangle of value onto the specified area. Define the 'real' dimensions of the big
% array, bottom-left corner of area to add. Must have dimensions of at
% least one cell in each direction in order to appear?
function [newArr] = addRect(baseArr, pos, smallDims, cellDims, val)
indeces = round(pos ./ cellDims)+1; % mfw matlab doesn't 0-index
lens = round(smallDims ./ cellDims);
tops = min(indeces + lens, size(baseArr));
bottoms = max(indeces, 0);
baseArr(bottoms(1):tops(1)-1, bottoms(2):tops(2)-1) = val * ones(tops-bottoms);
newArr = baseArr;
end