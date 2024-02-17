clc;
clear;

grain = 100;
pnum = 750;


% units should be kept in meter, second, kilogram, coloumb/volt, etc

% also uhh each particle will represent one micromole of ions?

fixedPots = zeros(grain+1, grain+1);
highPots = 15000;
lowPots = -400;
solids = [0.013, 0.00, 40*10^-6, highPots;...
    0.013, 0.0005, 40*10^-6, highPots;...
    0.015, 0.00, 40*10^-6, lowPots;...
    0.015, 0.0005, 40*10^-6, lowPots]; % x, y, mass, Volt in this case

%particles; % x, y, vx, vy, m, c
potentials = zeros(grain + 1,grain + 1); % evenly spaced in a grid
efieldx = zeros(grain + 1,grain + 1);
efieldy = zeros(grain + 1,grain + 1);
tstep = 0.001; % TODO make this adaptive, tstep <= 0.33 * spacing / mag(vel)

xmin = 0;
ymin = 0;
xmax = 0.05;
ymax = 0.02;
spacing = (xmax - xmin) / grain;



elecConstant = 1 / (4 * pi * 8.854*10^-12);

% pos_0 = init, pos_i = pos_i-1 + v_i-0.5 * dt
% v_0 = init, v_0.5 = v_0 + a_0 * dt / 2, v_i+0.5 = v_i-0.5 + a_i * dt
% ok so find efield distribution, then use that to find v_0.5
% so init is i = 0, first loop will be i = 1. so during first loop, i=1,
% pos will start with i=0 finding i=1, velocities will be 0.5 finding 1.5

solidDims = size(solids);

% for index = 1:solidDims(1)
%     fixedPots(round(solids(index, 1) / (xmax - xmin) * grain + 1), round(solids(index, 2) / (ymax - ymin) * grain + 1)) = solids(index, 4);
% end

fixedPots(14,1:20) = ones(1,20) * highPots;
fixedPots(44,1:20) = ones(1,20) * lowPots;

fixedPots(1,:) = 800 * ones(1,grain+1);
fixedPots(grain + 1,:) = 1 * ones(1,grain + 1);

spaceCharge = zeros(grain+1, grain+1);
oldSpace = spaceCharge;
newSpace = ones(size(spaceCharge));

spaceIters = 0;
maxSPaceIters = 100;

originalParticles = [ones(pnum, 1) * 0, rand(pnum, 1) * (ymax - ymin) + ymin, ones(pnum, 4) .* [50, 0, 39.948*10^-0, 9.65*10^-1]];

while max(max(abs(oldSpace - newSpace))) > 10^-20 && spaceIters < maxSPaceIters

    particles = originalParticles;

    oldSpace = newSpace;
    newSpace = zeros(size(oldSpace));
    
    oldPots = potentials;
    newPots = fixedPots;
    
    iters = 0;
    maxIters = 10000;
    while max(max(abs(newPots - oldPots))) > 0.0001 && iters < maxIters % check whether this tolerance is good, also make some iteration terminator
        oldPots = newPots;
        newPots = zeros(size(fixedPots));
        for xIndex = 1:(grain + 1)
            for yIndex = 1:(grain + 1)
                if fixedPots(xIndex, yIndex) == 0
                    xUp = xIndex + 1;
                    xDown = xIndex - 1;
                    yUp = yIndex + 1;
                    yDown = yIndex - 1;
                    % xDiff = 0;
                    % yDiff = 0;
                    % divisor = 0;
                    xDelta = ((xmax - xmin) / grain).^2;
                    yDelta = ((ymax - ymin) / grain).^2;
                    if xDown ~= 0 && xUp ~= grain + 2
                        
                    elseif xDown == 0
                        % xDiff = (-2 * oldPots(xUp, yIndex) + oldPots(xUp + 1, yIndex)) / xDelta;
                        % divisor = divisor + 1 / xDelta;
                        xDown = grain + 1;
                    else
                        % xDiff = (-2 * oldPots(xDown, yIndex) + oldPots(xDown - 1, yIndex)) / xDelta;
                        % divisor = divisor + 1 / xDelta;
                        xUp = 1;
                    end
                    if yDown ~= 0 && yUp ~= grain + 2
                        
                        % divisor = divisor - 2 / yDelta;
                    elseif yDown == 0
                        % yDiff = (-2 * oldPots(xIndex, yUp) + oldPots(xIndex, yUp + 1)) / yDelta;
                        % divisor = divisor + 1 / yDelta;
                        yDown = grain + 1;
                    else
                        % yDiff = (-2 * oldPots(xIndex, yDown) + oldPots(xIndex, yDown - 1)) / yDelta;
                        % divisor = divisor + 1 / yDelta;
                        yUp = 1;
                    end
                    xDiff = (oldPots(xDown, yIndex) + oldPots(xUp, yIndex)) / xDelta;
                    yDiff = (oldPots(xIndex, yDown) + oldPots(xIndex, yUp)) / yDelta;
                    divisor = -2 / xDelta - 2 / yDelta;
                    newPots(xIndex, yIndex) = (spaceCharge(xIndex, yIndex) - xDiff - yDiff) / divisor;
                end
            end
        end
        newPots = newPots + fixedPots;
        iters = iters + 1;
    end
    potentials = newPots;
    
    xDelta = (xmax - xmin) / grain;
    yDelta = (ymax - ymin) / grain;
    
    for xIndex = 1:(grain + 1)
        for yIndex = 1:(grain + 1)
            posx = (xIndex - 1) * (xmax - xmin) / (grain);
            posy = (yIndex - 1) * (ymax - ymin) / (grain);
            for solidIndex = 1:solidDims(1)
                % eAbs = elecConstant * solids(solidIndex, 4) / ((solids(solidIndex, 1) - posx).^2 + (solids(solidIndex, 2) - posy).^2);
                % angle = atan2(solids(solidIndex, 2) - posy, solids(solidIndex, 1) - posx);
                % efieldx(xIndex, yIndex) = efieldx(xIndex, yIndex) - eAbs * cos(angle);
                % efieldy(xIndex, yIndex) = efieldy(xIndex, yIndex) - eAbs * sin(angle);
                efieldx(xIndex, yIndex) = -(-potentials(mod(xIndex + grain - 1, grain + 1) + 1, yIndex) + potentials(mod(xIndex, grain + 1) + 1, yIndex)) / xDelta;
                efieldy(xIndex, yIndex) = -(-potentials(xIndex, mod(yIndex + grain - 1, grain + 1) + 1) + potentials(xIndex, mod(yIndex, grain + 1) + 1)) / yDelta;
            end
        end
    end
    
    pDims = size(particles);
    
    
    
    for pIndex = 1:pDims(1)
        xSubdiv = particles(pIndex,1) / (xmax - xmin) * grain;
        xFraction = xSubdiv - floor(xSubdiv); % kinda just abuse unit square properties
        ySubdiv = particles(pIndex,2) / (ymax - ymin) * grain;
        yFraction = ySubdiv - floor(ySubdiv); % add 1 to indeces bc FUCK matlab!
        accelX = particles(pIndex, 6) / particles(pIndex, 5) * lerp2(xSubdiv, xFraction, ySubdiv, yFraction, efieldx);
        accelY = particles(pIndex, 6) / particles(pIndex, 5) * lerp2(xSubdiv, xFraction, ySubdiv, yFraction, efieldy);
        particles(pIndex, 3) = particles(pIndex, 3) + accelX * tstep / 2;
        particles(pIndex, 4) = particles(pIndex, 4) + accelY * tstep / 2;
    end
    
    maxSteps = 5000;
    step = 0;
    time = 0;
    
    xTraces = 0;
    yTraces = 0;
    
    for pIndex = 1:pDims(1)
        step = 0;
        time = 0;
        temptraces = zeros(maxSteps, 3);
        temptraces(1,:) = [0, particles(pIndex, 1), particles(pIndex, 2)];
        while step < maxSteps && particles(pIndex,1) >= xmin && particles(pIndex,1) < xmax && particles(pIndex,2) >= ymin && particles(pIndex,2) < ymax
            xSubdiv = particles(pIndex,1) / (xmax - xmin) * grain;
            xFraction = xSubdiv - floor(xSubdiv); % kinda just abuse unit square properties
            ySubdiv = particles(pIndex,2) / (ymax - ymin) * grain;
            yFraction = ySubdiv - floor(ySubdiv); % add 1 to indeces bc FUCK matlab!

            
            accelX = particles(pIndex, 6) / particles(pIndex, 5) * lerp2(xSubdiv, xFraction, ySubdiv, yFraction, efieldx);
            accelY = particles(pIndex, 6) / particles(pIndex, 5) * lerp2(xSubdiv, xFraction, ySubdiv, yFraction, efieldy);

            xSubdiv = floor(xSubdiv);
            ySubdiv = floor(ySubdiv);

            newSpace(xSubdiv + 1, ySubdiv + 1) = newSpace(xSubdiv + 1, ySubdiv + 1) + particles(pIndex, 6) * (1 - xFraction) * (1 - yFraction) / elecConstant;
            newSpace(xSubdiv + 2, ySubdiv + 1) = newSpace(xSubdiv + 2, ySubdiv + 1) + particles(pIndex, 6) * (xFraction) * (1 - yFraction) / elecConstant;
            newSpace(xSubdiv + 1, ySubdiv + 2) = newSpace(xSubdiv + 1, ySubdiv + 2) + particles(pIndex, 6) * (1 - xFraction) * (yFraction) / elecConstant;
            newSpace(xSubdiv + 2, ySubdiv + 2) = newSpace(xSubdiv + 2, ySubdiv + 2) + particles(pIndex, 6) * (xFraction) * (yFraction) / elecConstant;

        
            tstep = min(0.001, spacing * getStepDivSpace(particles(pIndex,:)));
        
            particles(pIndex, 1) = particles(pIndex, 1) + particles(pIndex, 3) * tstep;
            particles(pIndex, 2) = particles(pIndex, 2) + particles(pIndex, 4) * tstep;
        
            particles(pIndex, 3) = particles(pIndex, 3) + accelX * tstep;
            particles(pIndex, 4) = particles(pIndex, 4) + accelY * tstep;
            
            step = step + 1;
            time = time + tstep;
            temptraces(step,:) = [time, particles(pIndex,1), particles(pIndex,2)];
        
        end
    
        while step <= maxSteps
            temptraces(step,:) = [time, particles(pIndex,1), particles(pIndex,2)];
            step = step + 1;
        end
    
        if(xTraces == 0)
            xTraces = temptraces(:,2);
            yTraces = temptraces(:,3);
        else
            xTraces = [xTraces, temptraces(:,2)];
            yTraces = [yTraces, temptraces(:,3)];
        end
    end
    
    G.tracefig = findobj(0,'name','sim');
    if isempty(G.tracefig)
        G.tracefig = figure();
    end
    figure(G.tracefig);
    clf;
    
    set(G.tracefig,...
        'Color', 'white',...
        'Menubar', 'figure',...
        'NumberTitle', 'off',...
        'Name', 'sim');
    colormap bone;
    
    set(axes(),...
        'Xlim', [xmin, xmax],...
        'Ylim', [ymin, ymax])
    
    plot(xTraces, yTraces)
    
    G.fieldfig = findobj(0,'name','field');
    if isempty(G.fieldfig)
        G.fieldfig = figure();
    end
    figure(G.fieldfig);
    clf;
    
    set(G.fieldfig,...
        'Color', 'white',...
        'Menubar', 'figure',...
        'NumberTitle', 'off',...
        'Name', 'field');
    colormap bone;
    
    % set(axes(),...
    %     'Xlim', [xmin, xmax],...
    %     'Ylim', [ymin, ymax])
    
    quiver(efieldx', efieldy', 1)
    
    
    G.potfig = findobj(0,'name','potentials');
    if isempty(G.potfig)
        G.potfig = figure();
    end
    figure(G.potfig);
    clf;
    
    set(G.potfig,...
        'Color', 'white',...
        'Menubar', 'figure',...
        'NumberTitle', 'off',...
        'Name', 'potentials');
    colormap bone;
    
    % set(axes(),...
    %     'Xlim', [xmin, xmax],...
    %     'Ylim', [ymin, ymax])
    
    pcolor(potentials')

    spaceCharges = newSpace;

    G.spacefig = findobj(0,'name','spacecharge');
    if isempty(G.spacefig)
        G.spacefig = figure();
    end
    figure(G.spacefig);
    clf;
    
    set(G.spacefig,...
        'Color', 'white',...
        'Menubar', 'figure',...
        'NumberTitle', 'off',...
        'Name', 'spacecharge');
    colormap bone;
    
    % set(axes(),...
    %     'Xlim', [xmin, xmax],...
    %     'Ylim', [ymin, ymax])
    
    pcolor(spaceCharges')

    
    fprintf('iteration %g\n', spaceIters)

    spaceIters = spaceIters + 1;

end

function [ret] = lerp2(x, xF, y, yF, arr)
bottom = arr(floor(x) + 1, floor(y) + 1) + (arr(floor(x) + 2, floor(y) + 1) - arr(floor(x) + 1, floor(y) + 1)) * xF;
top = arr(floor(x) + 1, floor(y) + 2) + (arr(floor(x) + 2, floor(y) + 2) - arr(floor(x) + 1, floor(y) + 2)) * xF;
ret = (bottom + (top - bottom) * yF);
end

function [stepsize] = getStepDivSpace(in)
stepsize = 1 / 3 / sqrt(in(3).^2 + in(4).^2);
end