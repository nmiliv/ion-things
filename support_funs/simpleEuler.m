function [initParticle] = simpleEuler(particle, tstep, forces)
accels = forces ./ particle(5);
particle(1:2) = particle(1:2) + particle(3:4) .* tstep;
particle(3:4) = particle(3:4) + accels .* tstep;
initParticle = particle;
end