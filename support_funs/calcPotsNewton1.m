% doesn't implement anyithing electron related


function [newPots, iters] = calcPotsNewton1(fixedPots, oldPots, spaceCharge, ...
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
diffs = ones(dims) * 2 * tol;

% fitness condition could be modified, currently runs until all values are
% moving slow enough, but could be run until all values are changing under
% a certain percentage or rate of convergence slows past a certain point.
while any(any(diffs > tol)) && iters < maxIters
    newPots = fixedPots;
    
    for xIndex = 1:dims(1)
        for yIndex = 1:dims(2)
            if fixedPots(xIndex, yIndex) == 0
                xDown = oldPots(max(1, xIndex - 1), yIndex); % mirror(xIndex - 1, 1, grain(1)+1);
                xUp = oldPots(min(dims(1), xIndex + 1), yIndex); % mirror(xIndex + 1, 1, grain(1)+1);
                yDown = oldPots(xIndex, mirror(yIndex - 1, 1, dims(2)));
                yUp = oldPots(xIndex, mirror(yIndex, 1, dims(2)));
                charge = spaceCharge(xIndex, yIndex);
                
                
                % pretty sure this should all be simple root
                newPots(xIndex, yIndex) = oldPots(xIndex, yIndex) - (charge + (xDown - 2*oldPots(xIndex, yIndex) + xUp)/cellDims(1).^2 + (yDown - 2*oldPots(xIndex, yIndex) + yUp)/cellDims(2).^2)/(-2*sum(cellDims.^-2));
                % also, seems to converge at about the same speed as Jacobi
                % at the moment, in terms of iterations. Real time
                % convergence is about 5x slower, maybe due to assigning
                % names to variables instead of direct lookup in oldPots?
            end
        end
    end
    diffs = abs(newPots - oldPots);
    oldPots = newPots;
    iters = iters + 1;
end

end