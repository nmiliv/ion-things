function [initParticle] = leapfrogInit1(particle, tstep, forces)
accels = forces ./ particle(5);
particle(3:4) = particle(3:4) + accels .* tstep/2;
particle(1:2) = particle(1:2);
particle(6) = tstep;
initParticle = particle;
end