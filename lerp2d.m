% takes (x, y) coords mapped to mesh coords and fractional distance between
% mesh coords. Returns the linear interpolation from the array
function [ret] = lerp2d(x, xF, y, yF, arr)
bottom = arr(floor(x) + 1, floor(y) + 1) + (arr(floor(x) + 2, floor(y) + 1) - arr(floor(x) + 1, floor(y) + 1)) * xF;
top = arr(floor(x) + 1, floor(y) + 2) + (arr(floor(x) + 2, floor(y) + 2) - arr(floor(x) + 1, floor(y) + 2)) * xF;
ret = (bottom + (top - bottom) * yF);
end