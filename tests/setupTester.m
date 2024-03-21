clc, clear;

grain = [200, 50];
maxes = [0.03, 0.01]; % max axial, radial distance
cellDeltas = maxes./grain;

elecTemp = 5; % electron volts
plasmaPot = 2266;
neutralPot = -0.1;
screenPot = 2241;
accelPot = -400;
sheathPot = plasmaPot - 0.5 * elecTemp;
currentBeamlet = 1*1e-4; % amps

fixedPots = zeros(grain+1);
% pretty values: plasma 1600, screen 2000, accel -100
% add fixed potentials here. a fixed potential of zero will indicate that
% there is no solid there and the potential should not be fixed (i.e. is a
% 'driven' potential instead of a 'driving' potential). So if you want to
% fix a potential to zero, just fix to some small number like 0.1.
% fixedPots = addRect(fixedPots, [0.005, 0.00], [0.002, 0.005], cellDeltas, 1);
fixedPots = addRect(fixedPots, [0.01, 0.00], [0.0004, 0.005], cellDeltas, screenPot);
fixedPots = addRect(fixedPots, [0.00, 0.00], [1.5*cellDeltas(1), maxes(2)+cellDeltas(2)], cellDeltas, plasmaPot);
fixedPots = addRect(fixedPots, [maxes(1)-0.001, 0.00], [0.001, maxes(2)+cellDeltas(2)], cellDeltas, neutralPot);
fixedPots = addRect(fixedPots, [0.0114, 0.00], [0.0008, 0.004], cellDeltas, accelPot);

maxPotIters = 50000; % potential convergence cutoff (50k is usually fine, but for meshes greater than 80k points will need to be even bigger)
potTolerance = 0.5 * 10^-3; % you get the point. 10^-1 usually decent

oldPots = fixedPots;

tic
[jacobiPots, jIters] = calcPotsJacobi2(fixedPots, oldPots, zeros(grain+1), plasmaPot, elecTemp, cellDeltas, potTolerance, maxPotIters);
toc
fprintf("jacobi, %g iters\n",jIters);
findPotErrors(jacobiPots, fixedPots, zeros(grain+1), cellDeltas);

tic
[newtonPots, nIters] = calcPotsNewton2(fixedPots, oldPots, zeros(grain+1), plasmaPot, elecTemp, cellDeltas, potTolerance, maxPotIters);
toc
fprintf("newton, %g iters\n",nIters);
findPotErrors(newtonPots, fixedPots, zeros(grain+1), cellDeltas);