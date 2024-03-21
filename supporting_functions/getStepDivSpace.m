% things can get innacurate if you travel more than one-third of a grid
% space in one time step. Returns the maximum allowable time step given the
% grid spacing.
function [stepsize] = getStepDivSpace(vx, vy, xjump, yjump)
stepsize = min(abs([xjump/3/vx, yjump/3/vy]));
end