initV = 1;
initR = 0.3;
simTime = 20;
maxSteps = 200;

particleLFOld = [1 - initR, 1, 0, initV, 1];
particleLF1 = [1 - initR, 1, 0, initV, 1, 0];
particleEuler = [1 - initR, 1, 0, initV, 1];
particleLF2 = [1 - initR, 1, 0, initV, 1];
maxes = [2,2];
grain = [15, 15];
cellDeltas = maxes./(grain - 1);
forcesX = zeros(grain);
forcesY = zeros(grain);
for indexX = 1:grain(1)
    for indexY = 1:grain(2)
        pos = [indexX - 1, indexY - 1] .* cellDeltas;
        diag = sqrt(sum((pos-maxes/2).^2));
        normPos = (pos - maxes/2)./diag;
        forcesX(indexX, indexY) = -1 * initV.^2/initR * normPos(1);
        forcesY(indexX, indexY) = -1 * initV.^2/initR * normPos(2);
    end
end

G.fig = findobj(0, 'name', 'figs');
if isempty(G.fig)
    G.fig = figure();
end
figure(G.fig);
clf;

set(G.fig,...
    'Color', 'white',...
    'Menubar', 'figure', ...
    'NumberTitle', 'off', ...
    'Name', 'figs');
tiledlayout(1, 2);
nexttile;
ax1 = quiver(forcesX', forcesY');
nexttile;




time = 0;
step = 0;
xtraces = zeros(maxSteps,8);
ytraces = zeros(maxSteps,8);
while time < simTime && step < maxSteps
    time = time + 0.1;
    step = step + 1;
    xtraces(step,1) = -initR*cos(initV/initR*time)+1;
    ytraces(step,1) = initR*sin(initV/initR*time)+1;
end

time = 0;
step = 0;
particle = particleLFOld;
[subdivs, fracs] = subFracs(particle(1:2), cellDeltas);
particle(1:4) = leapfrogInitOld(particle, getStepDivSpace(particle(3),particle(4),cellDeltas(1),cellDeltas(2)), [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesX), lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesY)]);
while time < simTime && step < maxSteps
    tstep = getStepDivSpace(particle(3), particle(4), cellDeltas(1), cellDeltas(2));
    [subdivs, fracs] = subFracs(particle(1:2), cellDeltas);
    particle(1:4) = leapfrogStepOld(particle, getStepDivSpace(particle(3),particle(4),cellDeltas(1),cellDeltas(2)), [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesX), lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesY)]);
    step = step + 1;
    time = time + tstep;
    xtraces(step,2) = particle(1);
    ytraces(step,2) = particle(2);
end

time = 0;
step = 0;
particle = particleLF1;
[subdivs, fracs] = subFracs(particle(1:2), cellDeltas);
particle = leapfrogInit1(particle, getStepDivSpace(particle(3),particle(4),cellDeltas(1),cellDeltas(2)), [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesX), lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesY)]);
while time < simTime && step < maxSteps
    tstep = getStepDivSpace(particle(3), particle(4), cellDeltas(1), cellDeltas(2));
    [subdivs, fracs] = subFracs(particle(1:2), cellDeltas);
    particle = leapfrogStep1(particle, getStepDivSpace(particle(3),particle(4),cellDeltas(1),cellDeltas(2)), [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesX), lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesY)]);
    step = step + 1;
    time = time + tstep;
    xtraces(step,3) = particle(1);
    ytraces(step,3) = particle(2);
end

time = 0;
step = 0;
particle = particleEuler;
while time < simTime && step < maxSteps
    tstep = getStepDivSpace(particle(3), particle(4), cellDeltas(1), cellDeltas(2));
    [subdivs, fracs] = subFracs(particle(1:2), cellDeltas);
    particle = simpleEuler(particle, getStepDivSpace(particle(3),particle(4),cellDeltas(1),cellDeltas(2)), [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesX), lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesY)]);
    step = step + 1;
    time = time + tstep;
    xtraces(step,4) = particle(1);
    ytraces(step,4) = particle(2);
end

% LF2 method
particle = [1 - initR, 1, 0, initV, 1];
time = 0;
step = 0;
while time < simTime && step < maxSteps
    tstep = getStepDivSpace(particle(3), particle(4), cellDeltas(1), cellDeltas(2));
    [subdivs, fracs] = subFracs(particle(1:2), cellDeltas);
    particle = leapfrogStep2(particle, getStepDivSpace(particle(3),particle(4),cellDeltas(1),cellDeltas(2)), [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesX), lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesY)]);
    step = step + 1;
    time = time + tstep;
    xtraces(step,5) = particle(1);
    ytraces(step,5) = particle(2);
end

% RK4 method
particle = [1 - initR, 1, 0, initV, 1];
time = 0;
step = 0;
while time < simTime && step < maxSteps
    tstep = getStepDivSpace(particle(3), particle(4), cellDeltas(1), cellDeltas(2));
    particle = stepRK4(particle, tstep, forcesX, forcesY, cellDeltas);
    step = step + 1;
    time = time + tstep;
    xtraces(step,6) = particle(1);
    ytraces(step,6) = particle(2);
end


% define potential energy as distance from (1,1), try energy conservation
% with some of the leapfrog methods

time = 0;
step = 0;
particle = particleLFOld;
initE = initR + 0.5*initV.^2;
[subdivs, fracs] = subFracs(particle(1:2), cellDeltas);
particle(1:4) = leapfrogInitOld(particle, getStepDivSpace(particle(3),particle(4),cellDeltas(1),cellDeltas(2)), [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesX), lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesY)]);
while time < simTime && step < maxSteps
    tstep = getStepDivSpace(particle(3), particle(4), cellDeltas(1), cellDeltas(2));
    [subdivs, fracs] = subFracs(particle(1:2), cellDeltas);
    particle(1:4) = leapfrogStepOld(particle, getStepDivSpace(particle(3),particle(4),cellDeltas(1),cellDeltas(2)), [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesX), lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesY)]);
    particle(3:4) = particle(3:4) .* sqrt(2 * (initE - sqrt(sum((particle(1:2) - 1).^2))) / sum(particle(3:4).^2));
    step = step + 1;
    time = time + tstep;
    xtraces(step,7) = particle(1);
    ytraces(step,7) = particle(2);
end

particle = [1 - initR, 1, 0, initV, 1];
time = 0;
step = 0;
initE = initR + 0.5*initV.^2;
while time < simTime && step < maxSteps
    tstep = getStepDivSpace(particle(3), particle(4), cellDeltas(1), cellDeltas(2));
    [subdivs, fracs] = subFracs(particle(1:2), cellDeltas);
    particle = leapfrogStep2(particle, getStepDivSpace(particle(3),particle(4),cellDeltas(1),cellDeltas(2)), [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesX), lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesY)]);
    particle(3:4) = particle(3:4) .* sqrt(2 * (initE - sqrt(sum((particle(1:2) - 1).^2))) / sum(particle(3:4).^2));
    step = step + 1;
    time = time + tstep;
    xtraces(step,8) = particle(1);
    ytraces(step,8) = particle(2);
end


ax2 = plot(xtraces, ytraces);
legend("true","LF Old", "LF 1", "Euler", "LF2", "RK4", "LF Old CoE", "LF2 CoE");


% stress testing:
stress = 10000;
fprintf("RK4: ");
tic
for index = 1:stress
    particle = [1 - initR, 1, 0, initV, 1];
    time = 0;
    step = 0;
    while step < maxSteps
        tstep = getStepDivSpace(particle(3), particle(4), cellDeltas(1), cellDeltas(2));
        particle = stepRK4(particle, tstep, forcesX, forcesY, cellDeltas);
        step = step + 1;
        time = time + tstep;
        xtraces(step,6) = particle(1);
        ytraces(step,6) = particle(2);
    end
end
toc

fprintf("LF old: ");
tic
for index = 1:stress
    time = 0;
    step = 0;
    particle = particleLFOld;
    initE = initR + 0.5*initV.^2;
    [subdivs, fracs] = subFracs(particle(1:2), cellDeltas);
    particle(1:4) = leapfrogInitOld(particle, getStepDivSpace(particle(3),particle(4),cellDeltas(1),cellDeltas(2)), [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesX), lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesY)]);
    while time < simTime && step < maxSteps
        tstep = getStepDivSpace(particle(3), particle(4), cellDeltas(1), cellDeltas(2));
        [subdivs, fracs] = subFracs(particle(1:2), cellDeltas);
        particle(1:4) = leapfrogStepOld(particle, getStepDivSpace(particle(3),particle(4),cellDeltas(1),cellDeltas(2)), [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesX), lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesY)]);
        particle(3:4) = particle(3:4) .* sqrt(2 * (initE - sqrt(sum((particle(1:2) - 1).^2))) / sum(particle(3:4).^2));
        step = step + 1;
        time = time + tstep;
        xtraces(step,7) = particle(1);
        ytraces(step,7) = particle(2);
    end
end
toc

fprintf("LF2: ");
tic
for index = 1:stress
    particle = [1 - initR, 1, 0, initV, 1];
    time = 0;
    step = 0;
    initE = initR + 0.5*initV.^2;
    while time < simTime && step < maxSteps
        tstep = getStepDivSpace(particle(3), particle(4), cellDeltas(1), cellDeltas(2));
        [subdivs, fracs] = subFracs(particle(1:2), cellDeltas);
        particle = leapfrogStep2(particle, getStepDivSpace(particle(3),particle(4),cellDeltas(1),cellDeltas(2)), [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesX), lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesY)]);
        particle(3:4) = particle(3:4) .* sqrt(2 * (initE - sqrt(sum((particle(1:2) - 1).^2))) / sum(particle(3:4).^2));
        step = step + 1;
        time = time + tstep;
        xtraces(step,8) = particle(1);
        ytraces(step,8) = particle(2);
    end
end
toc