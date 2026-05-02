function y = Trap(t, G, tStart, tDuration, tRampTime)
%TRAP Summary of this function goes here
%   Detailed explanation goes here

 y = RampU(t,G,tStart,tRampTime) + SquarePulse(t,G,tStart + tRampTime,tDuration-tRampTime) + RampD(t,G,tStart+tDuration,tRampTime);

end
