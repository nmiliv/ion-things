function [subdivs, fracs] = subFracs(posXY, cellDims)
subdivs = posXY ./ cellDims;
fracs = subdivs - floor(subdivs);
subdivs = floor(subdivs);
end