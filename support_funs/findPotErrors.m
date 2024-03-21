function [errors] = findPotErrors(potentials, fixedPots, spacecharge, cellDims)
dims = size(potentials);
errors = zeros(dims);
divisor = -2 * sum(cellDims.^-2);
for xIndex = 1:dims(1)
    for yIndex = 1:dims(2)
        if fixedPots(xIndex, yIndex) == 0
            xDown = max(1, xIndex - 1); % mirror(xIndex - 1, 1, grain(1)+1);
            xUp = min(dims(1), xIndex + 1); % mirror(xIndex + 1, 1, grain(1)+1);
            yDown = mirror(yIndex - 1, 1, dims(2));
            yUp = mirror(yIndex, 1, dims(2));
            charge = spacecharge(xIndex, yIndex);
            %{
            errors(xIndex, yIndex) = abs((potentials(xDown, yIndex) - 2 * potentials(xIndex, yIndex) + potentials(xUp, yIndex)) / cellDims(1).^2 + ...
                (potentials(xIndex, yDown) - 2 * potentials(xIndex, yIndex) + potentials(xIndex, yUp)) / cellDims(2).^2 + ...
                charge);
            %}
            errors(xIndex, yIndex) = abs((-charge...
                - (potentials(xDown, yIndex) + potentials(xUp, yIndex))/(cellDims(1)^2)...
                - (potentials(xIndex, yDown) + potentials(xIndex, yUp))/(cellDims(2)^2))/divisor...
                - potentials(xIndex, yIndex));

        end
    end
end
fprintf("Maximum error: %g\n", max(max(errors)));
end