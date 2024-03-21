function [initParticle] = leapfrogStep2(particle, tstep, forces)
accels = forces ./ particle(5);
particle(1:2) = particle(1:2) + particle(3:4) .* tstep + accels .* tstep.^2 / 2;
particle(3:4) = particle(3:4) + accels .* tstep;
initParticle = particle;
end