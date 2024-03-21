function [initParticle] = leapfrogStepOld(particle, tstep, forces)
accels = forces ./ particle(5);
initParticle(1:2) = particle(1:2) + particle(3:4) .* tstep; % since velocities are for n = t+1/2, position n = t+1 gets updated using 'old' velocity
initParticle(3:4) = particle(3:4) + accels .* tstep;
end