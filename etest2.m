clc;
clear;

% initialization paramters
% important to note that x->rows, y->cols. so all arrs must be transposed
% to 'display' properly. Also means that (x,y) translates nicely to
% (row,col) format.
pnum = 500; % number of particles each iteration
grain = [500, 100]; % (x,y) ordered pair of cells in each direction
% we'll assume that eveything starts at the origin, quadrant 1.
maxes = [0.02, 0.002]; % max axial, radial distance
cellDeltas = maxes./grain;
maxSimSteps = 5000; % TODO wait this should be reactive to mesh size bc we can prove a particle mut take at least 3*n steps to cross the sim
% this'll need to be bigger than 3*n tho, bc particles might loop around or
% something idk. honestly this hard limit should get removed.


elecTemp = 5; % electron volts
plasmaPot = 2266;
neutralPot = -0.1;
screenPot = 2241;
accelPot = -400;
sheathPot = plasmaPot - 0.5 * elecTemp;
currentBeamlet = 0.5*1e-4; % amps

fixedPots = zeros(grain+1);
% pretty values: plasma 1600, screen 2000, accel -100
% add fixed potentials here. a fixed potential of zero will indicate that
% there is no solid there and the potential should not be fixed (i.e. is a
% 'driven' potential instead of a 'driving' potential). So if you want to
% fix a potential to zero, just fix to some small number like 0.1.
% fixedPots = addRect(fixedPots, [0.005, 0.00], [0.002, 0.005], cellDeltas, 1);
fixedPots = addRect(fixedPots, [0.005, 0.00], [0.0004, 0.001], cellDeltas, screenPot);
fixedPots = addRect(fixedPots, [0.00, 0.00], [1*cellDeltas(1), maxes(2)+cellDeltas(2)], cellDeltas, plasmaPot);
fixedPots = addRect(fixedPots, [maxes(1)-0.001, 0.00], [0.001, maxes(2)+cellDeltas(2)], cellDeltas, neutralPot);
fixedPots = addRect(fixedPots, [0.0064, 0.00], [0.0008, 0.0014], cellDeltas, accelPot);


%{
so I guess it's time to talk about how the nodes work bc that's actually
important now. Pcolor is a bit misleading here, since it colors a cell
based on the bottom-left corner and ignore the topmost/leftmost cells.
Don't trust pcolor.

Here's some nodes, marked as X:
X--X--X
|  |  |
|  |  |
X--X--X
|  |  |
|  |  |
X--X--X

If a node has a fixedPot, that mean's it's probably a solid (make sure to
turn off solid checking in known plasma/space regions -- a bit of a hack
but honestly not too bad imo??) The solid region for the center node is
highlighted with o's:

X--X--X
|  |  |
| ooo |
X-oXo-X
| ooo |
|  |  |
X--X--X

basically, the node is at the center of a cell rather than the nodes
defining the corners of the cells (node-is-cell vs. node-defines-cell). Is
this a little innacurate? Yeah probably. But as long as the mesh is small
enough it should be close enough lmao.

TODO shit check if the bohm velocity is respected in the half-step velocity
init thing
%}



% actually process initialization
% [x, y, vx, vy, mass/charge]
massCharge = (2.1802*10^-6)/(1.602*10^-0); % mass / charge
vBohm = 1 * sqrt(elecTemp / massCharge);
ogParticles = [ones(pnum, 1) * 2 * cellDeltas(1), ...
    linspace(0,0.9999999999*maxes(2),pnum)', ...
    ones(pnum, 3) .* [vBohm, 0, massCharge]];
% also it's supposed to be m = 6.6335*10^-10, c = 1.602*10^-3 for argon but uhh in
% the name of pretty graphics we're gonna ignore that
% pretty results: m= 4*10^1, c = 9.6*10^-2, v = 5
% argon m = 6.6335*10^-7, c = 1.602*10^0
% initial velocity should be bohm velocity, which should probably be
% actually calculated instead of guesstimated.
% v_bohm = f_bohm * sqrt( q*T_e / m_i )
% f can be taken as 1, q is electric charge, T is electron temperature in
% electron volts, m is the mass of the ion

% want to make sure we all always start with the same particles so that
% changes in space charge are due to particle paths changing, not initial
% conditions changing
% TODO make sure that these values actually make sense, implement double
% ions, etc.


% this isn't going to worry about grid erosion, so defining an array of
% solids or driving grid params for that won't happen yet
maxSpaceIters = 1000; % maximum iterations for calculating space charge, realistically would take several hours or days of sim time to reach this
spaceTolerance = 10^2; % maximum allowable change between space charge iterations, this is typically on the order of 10^5
maxPotIters = 50000; % potential convergence cutoff (50k is usually fine, but for meshes greater than 80k points will need to be even bigger)
potTolerance = 0.5 * 10^-1; % you get the point. 10^-1 usually decent


% simulation paramters
spaceIters = 0;
currentParticle = currentBeamlet / pnum; % TODO this assumes all particles are singly charged!

% physical constants
eps0 = 8.854*10^-12; % epsilon naught
elecConstant = 1 / 4 / pi / eps0;
boltzmann = 1.38065 * 10^-23;

lenDebye = eps0 * boltzmann * elecTemp / (1.92*10^-19)^2;

% arrays that will be used to keep track of things later
newSpace = zeros(grain+1);
oldSpace = (newSpace + 1) * 2 * spaceTolerance; % just need to make sure this isn't the same as newspace, will get reset later (i.e. tricking the while loop into running the first time)
oldPots = fixedPots;
newPots = zeros(grain+1);
efieldx = zeros(grain+1);
efieldy = zeros(grain+1);


xTraces = zeros(maxSimSteps, pnum);
yTraces = zeros(maxSimSteps, pnum);

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

% TODO there should be some way of calculating the sheath and
% neutralization surface. also maybe dynamically changing the size of the
% simulation to account for this stuff

% main process loop
while max(max(abs(oldSpace - newSpace))) > spaceTolerance && spaceIters < maxSpaceIters % TODO make space tolerance adaptive
    
    % itinerary: use space charges and fixed voltages to find potential
    % field, turn that into e-field, track particles through that to update
    % the space charge
    particles = ogParticles;

    oldSpace = newSpace;
    % particles will behave based on oldspace, particle behaviour will be used to inform newspace
    newSpace = zeros(grain+1);

    
    [newPots, potIters] = calcPotsJacobi2(fixedPots, oldPots, oldSpace, plasmaPot, elecTemp, cellDeltas, potTolerance, maxPotIters);

    fprintf('Potential iterations: %g\n', potIters)

    oldPots = newPots;

    % calculate efield
    % ugh this could probably be parallelized in some way as well but
    % whatever
    for xIndex = 1:grain(1)+1
        for yIndex = 1:grain(2)+1
            % TODO yeah i think a hard clip would be better than a mirror
            % here
            xDown = max(xIndex(1) - 1, 1);
            xUp = min(xIndex(1) + 1, grain(1) + 1);
            yDown = mirror(yIndex - 1, 1, grain(2) + 1);
            yUp = mirror(yIndex + 1, 1, grain(2) + 1);
            efieldx(xIndex, yIndex) = -(-oldPots(xDown, yIndex) + oldPots(xUp, yIndex)) / cellDeltas(1);
            efieldy(xIndex, yIndex) = -(-oldPots(xIndex, yDown) + oldPots(xIndex, yUp)) / cellDeltas(2);
            
        end
    end

    % track particles
    % TODO ok this can def be parallelized
    % but imma just copy the old code for now
    % TODO ok actually look into how multithreading in matlab works aaaaa
    % (tic and toc, which matrix operations are default multithreaded, how
    % to force a multithread, etc)
    for pIndex = 1:pnum
        step = 0;
        time = 0;
        tempTraces = zeros(maxSimSteps, 3); % ok wait we could totally get rid of maxSimSteps if we don't plot the individual traces of the particles
        % or parelellizing will also solve the maxSimSteps issue
        particle = particles(pIndex,:);
        temptraces(1,:) = [particle(1:2)];
        

        % currently the particles are advanced using the 'leapfrog' method,
        % which is a slight improvement over the euler method while adding
        % basically no complexity. Essentially the velocity between time
        % steps is used (i.e. average velocity) instead of velocity at some
        % specific time step. this means the velocity values have to be
        % 'primed' by setting them to the n=1/2 value
        tstep = getStepDivSpace(particle(3), particle(4), cellDeltas(1), cellDeltas(2));
        [subdivs, fracs] = subFracs(particle(1:2), cellDeltas);
        accels = [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), efieldx),...
            lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), efieldy)] ./ particle(5);
        % my ass should NOT be out debugging this late T.T
        % there's allegedly some indexing issue out here
        % update from a couple weeks later: I think the indexing issue is
        % fixed. Honestly can't remember what that was about.
        particle(3:4) = particle(3:4) + accels * tstep / 2;
        % TODO still not sure how half-stepping works with the whole
        % adaptive tstep thing.

        % fprintf('Running particle %g\n', pIndex)
        while (step < maxSimSteps) && all(particle(1:2) >= 0) && all(particle(1:2) < maxes)
            [subdivs, fractions] = subFracs(particle(1:2), cellDeltas);
            accels = [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), efieldx),...
                lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), efieldy)] ./ particle(5);

            % TODO check for collisions here!
            if subdivs(1) > 1 && subdivs(1) < grain(1) - 2 && ...
                    fixedPots(subdivs(1) + 1 + round(fracs(1)), subdivs(2) + 1 + round(fracs(2))) ~= 0
                % fprintf("collision!\n");
                break;

            end
            
            % honestly haven't the faintest clue why this is getting done
            % in this order. also there's probably some fancier matrix way
            % of doing this operation but i have a morning class tomorrow
            vCell = cellDeltas(1) * cellDeltas(2);
            newSpace(subdivs(1) + 1, subdivs(2) + 1) = newSpace(subdivs(1) + 1, subdivs(2) + 1) + currentParticle * tstep / vCell * (1 - fractions(1)) * (1 - fractions(2)) / eps0;
            newSpace(subdivs(1) + 2, subdivs(2) + 1) = newSpace(subdivs(1) + 2, subdivs(2) + 1) + currentParticle * tstep / vCell * (fractions(1)) * (1 - fractions(2)) / eps0;
            newSpace(subdivs(1) + 1, subdivs(2) + 2) = newSpace(subdivs(1) + 1, subdivs(2) + 2) + currentParticle * tstep / vCell * (1 - fractions(1)) * (fractions(2)) / eps0;
            newSpace(subdivs(1) + 2, subdivs(2) + 2) = newSpace(subdivs(1) + 2, subdivs(2) + 2) + currentParticle * tstep / vCell * (fractions(1)) * (fractions(2)) / eps0;
            % so I'm not sure if this should be defined as the total
            % current entering the sheath or just the current that makes it
            % out through the beamlet

            tstep = getStepDivSpace(particle(3), particle(4), cellDeltas(1), cellDeltas(2));
            
            particle(1:2) = particle(1:2) + particle(3:4) * tstep;
            mirrorPos = mirror(particle(2), 0, maxes(2));
            particle(3:4) = particle(3:4) + accels * tstep;
            if mirrorPos ~= particle(2)
                particle(2) = mirrorPos;
                particle(4) = -particle(4);
            end

            
            
            
            % ugh things are still broken here
            % TODO add energy conservation step!
            % ok so find maximum time step, also need to init the velocities of the particles somewhere, also make sure we're conserving energy
            % anyway it's like 1am and I need to pass the hell out
            % pretty sure i finished most of that stuff?

            step = step + 1;
            time = time + tstep;
            temptraces(step,:) = [particle(1:2)];

        end

        while step <= maxSimSteps
            temptraces(step,:) = [particle(1:2)];
            step = step + 1;
        end
        xTraces(:,pIndex) = temptraces(:,1);
        yTraces(:,pIndex) = temptraces(:,2);
    end

    % draw figures
    % ugh the subplot thing....
    figure(G.fig);
    clf;
    
    tiledlayout(2,2);
    nexttile
    plot(xTraces, yTraces)
    nexttile
    quiver(efieldx', efieldy')
    nexttile
    pcolor(oldPots')
    nexttile
    pcolor(newSpace')

    fprintf('Space iteration: %g\n',spaceIters)
    spaceIters = spaceIters + 1;

end