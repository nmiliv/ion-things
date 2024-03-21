function [initParticle] = leapfrogInitOld(particle, tstep, forces)
accels = forces ./ particle(5);
initParticle(3:4) = particle(3:4) + accels .* tstep/2;
initParticle(1:2) = particle(1:2);
end