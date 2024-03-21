% similar to wrap, but mirrors the input along the bound if it is out of
% bounds
function [mirror] = mirror(in, inmin, inmax)
clamp = max(inmin, min(inmax, in)); % clamp the value to the given range
bounce = in - clamp; % how far off we were from the range
mirror = clamp - bounce; % bounce back on the clamp!
end