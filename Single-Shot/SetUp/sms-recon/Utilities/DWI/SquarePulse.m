function y = SquarePulse(t,G,tStart,tDuration)
%SQUAREPULSE Summary of this function goes here
%   Detailed explanation goes here

sizesT = length(t);
sizesG = length(G); 
y = repmat(G,1,sizesT) .* repmat(Pulse(t - tStart) - Pulse(t - (tStart + tDuration)),sizesG,1);

end

