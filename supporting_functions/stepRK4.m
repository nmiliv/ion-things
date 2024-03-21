function [initParticle] = stepRK4(particle, tstep, forcesX, forcesY, cellDeltas) % inputs might be a little redundant but throwing the kitchen sink at this for now
% holy cow, can't believe i'm actually writing an RK4 algorithm

% oh hmm RK4 is for a first order DE, and we've got a second order one
% aight so classic RK4 would be for a velocity field, take a particle,
% take an euler half step, take another half step (but starting at original
% start point), and then take a full step (starting from original start
% point). Then take a weighted average of all the steps.

% what a second order thing might look like? Take an euler half step (get
% new velocity, move to a new position). Take another euler half step...
% oh hmm checked stack overflow lemme rethink this
% F(x) = ma --> x'' = F(x)/m --> x'' - F(x)/m = 0
% proxy u, x' = u and u' = F(x)/m
% initial conditions: x is whatever pos we got, u is velocity
% k and l used as RK4 constants

% F(x)/m is sorta like the 'source' function here, i think?

[subdivs, fracs] = subFracs(particle(1:2), cellDeltas);
l1 = tstep * particle(3:4);
k1 = tstep * [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesX), lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesY)]./particle(5);

[subdivs, fracs] = subFracs(particle(1:2) + 0.5*l1, cellDeltas);
l2 = tstep * (particle(3:4) + 0.5 * k1);
k2 = tstep * [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesX), lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesY)]./particle(5);

[subdivs, fracs] = subFracs(particle(1:2) + 0.5*l2, cellDeltas);
l3 = tstep * (particle(3:4) + 0.5 * k2);
k3 = tstep * [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesX), lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesY)]./particle(5);

[subdivs, fracs] = subFracs(particle(1:2) + l3, cellDeltas);
l4 = tstep * (particle(3:4) + k3);
k4 = tstep * [lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesX), lerp2d(subdivs(1), fracs(1), subdivs(2), fracs(2), forcesY)]./particle(5);

particle(1:4) = particle(1:4) + ([l1,k1] + 2.*[l2,k2] + 2.*[l3,k3] + [l4,k4])./6;
initParticle = particle;
end