% overall would not recommend newton method in most cases. Convergence is
% mostly limited by communication time across the mesh so extra
% computational cost of newton doesn't make sense (convergence is
% iteration-limited).

function [newPots, iters] = calcPotsNewton2(fixedPots, oldPots, spaceCharge, ...
    upstreamPot, elecTemp, cellDims, tol, maxIters)
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
                differential = -2*sum(cellDims.^-2);
                if oldPots(xIndex, yIndex) > upstreamPot
                    charge = charge - charge*(1 + (oldPots(xIndex, yIndex) - upstreamPot)/elecTemp);
                    differential = differential - charge/elecTemp;
                else
                    differential = differential - charge*exp((oldPots(xIndex, yIndex) - upstreamPot)/elecTemp)/elecTemp;
                    charge = charge - charge*exp((oldPots(xIndex, yIndex) - upstreamPot)/elecTemp);
                end
                
                % pretty sure this should all be simple root
                newPots(xIndex, yIndex) = oldPots(xIndex, yIndex) - (charge + (xDown - 2*oldPots(xIndex, yIndex) + xUp)/cellDims(1).^2 + (yDown - 2*oldPots(xIndex, yIndex) + yUp)/cellDims(2).^2)/(differential);
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