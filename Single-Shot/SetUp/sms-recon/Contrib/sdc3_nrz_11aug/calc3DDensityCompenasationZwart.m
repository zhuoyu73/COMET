function [dcf] = calc3DDensityCompenasation(kx,ky,kz)
%calc3DDensityCompensation Calculate 3D compensation function for an
%arbitrary trajectory. Uses sdc3_MAT.c and sdc3grid_kernel.c from Jim
%Pipe's group and Nick Zwart.
%   Inputs: 
%       kx - column vector of kx trajectory [-N/2:N/2-1]
%       ky - column vector of ky trajectory [-N/2:N/2-1]
%       kz - column vector of kz trajectory [-Nz/2:Nz/2-1]
%
%   Output:
%       dcf - density compensation function


% Settings from testmex.m
numIter = 25;
osf     = 2.1;
verbose = 1;

% Prepping inputs for input to sdc3_MAT. Important to avoid segfaulting

maxKx = max(kx);
minKx = min(kx);
maxKy = max(ky);
minKy = min(ky);
maxKz = max(kz);
minKz = min(kz);

maxes = [maxKx,maxKy,maxKz,abs(minKx),abs(minKy),abs(minKz)];
effMtx = 2*max(maxes);
effMtxz = 2*max([maxKz,abs(minKz)]);

crds = cat(ndims(kx)+1,kx./effMtx,ky./effMtx,kz./effMtxz);
crds = permute(crds,[ndims(crds),1:ndims(crds)-1]);
% Note that the array that goes into sdc3_MAT must be [3 x Npts]! Or you
% wil segfault and kill MATLAB. BEEEEWAAAREEE!!!
dcf = sdc3_MAT(crds,numIter,effMtx,verbose,osf);



end

