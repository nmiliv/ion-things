% doesn't implement anyithing electron related

% also this seems to have magically improved on the old version while also
% decreasing convergence time. It's literally the same code, but the values
% it's returning seem to fit a bit better. Def still not perfect though,
% but less sensitive to initial conditions for free potentials
function [newPots, iters] = calcPotsJacobi1(fixedPots, oldPots, spaceCharge, ...
    cellDims, tol, maxIters)
dims = size(fixedPots);
if any(dims ~= size(oldPots)) || ...
    any(dims ~= size(spaceCharge))
    error('input arrays should be of equal size');
end

% if any(any(oldPots ~= fixedPots))
%     fprintf("calcPotsJacobi1: recevied unique oldPots\n")
% end
iters = 0;

newPots = fixedPots;
divisor = -2 * sum(cellDims.^-2);
diffs = ones(dims) * 2 * tol;

% fitness condition could be modified, currently runs until all values are
% moving slow enough, but could be run until all values are changing under
% a certain percentage or rate of convergence slows past a certain point.
while any(any(diffs > tol)) && iters < maxIters
    newPots = fixedPots;
    
    for xIndex = 1:dims(1)
        for yIndex = 1:dims(2)
            if fixedPots(xIndex, yIndex) == 0
                xDown = max(1, xIndex - 1); % mirror(xIndex - 1, 1, grain(1)+1);
                xUp = min(dims(1), xIndex + 1); % mirror(xIndex + 1, 1, grain(1)+1);
                yDown = mirror(yIndex - 1, 1, dims(2));
                yUp = mirror(yIndex, 1, dims(2));
                charge = spaceCharge(xIndex, yIndex);
                    
                xDiff = (oldPots(xDown, yIndex) + oldPots(xUp, yIndex)) / cellDims(1).^2;
                yDiff = (oldPots(xIndex, yDown) + oldPots(xIndex, yUp)) / cellDims(2).^2;
                newPots(xIndex, yIndex) = (- charge - xDiff - yDiff) / divisor;
            end
        end
    end
    diffs = abs(newPots - oldPots);
    oldPots = newPots;
    iters = iters + 1;
end

end