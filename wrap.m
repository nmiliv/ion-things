% if the input is within bounds, returns the input. if it is out of bounds,
% wraps it around to the other side of allowable bounds. Might only be
% defined for positive bounds and input less than one range out of bounds
function [wrap] = wrap(in, min, max)
wrap = mod(in + max - 2 * min + 1, max - min + 1) + min;
% basically transform such that min is at 0, max is at max-min+1 (range).
% +1 is for inclusivity (republicans hate this one simple trick!)* then
% modulo using the range, followed by transform back to the original range.
% 
% *this joke is about pearl clutching. I am not saying all republicans will
% pearl clutch whenever they see the word "inclusive," even in a
% mathematical context. I have yet to see a political party with no pearl
% clutchers, so neither am I saying this is an issue exclusive to the
% republican party.
end