function [kx,ky] = genvdkspace(FOV,N,ld,nint,gamp,gslew,tsamp,rotamount,rev_flag,gts);
%function kx,ky] = genvdkspace(FOV,N,ld,nint,gamp,gslew,tsamp,rotamount,rev_flag,gts);
% This function will generate the proper length of k-space
% trajectory.  It linearly interpolates the output of
% genspiral to the correct length and takes care of the 
% rotations for the interleaves.
% ld is the length of the data
% nint is the number of interleaves
% Brad Sutton, Bioengineering, UIUC
flag = 0;   %auto determine number of k-space points
            % just input ld = 0.3

if ~(exist('rotamount','var'))
  rotamount = 0;
end
if ~(exist('rev_flag','var'))
  rev_flag = 0;
end
if ~(exist('gts','var'))
  gts = 10e-6;
end


nk = ld/nint;
if round(nk) ~= nk
sprintf('Input should have num data pts/number of interleaves must be int')
end

if (nk == 0)
     flag = 1;
end


%dts = 4e-6;    %5e-6
[kxi,kyi] = genspivd_Kim(FOV,N,nint,gamp,gslew,gts);


kxt=interp1([0:gts:gts*length(kxi)-gts],kxi,[0:tsamp:gts*length(kxi)-tsamp])';
kyt=interp1([0:gts:gts*length(kyi)-gts],kyi,[0:tsamp:gts*length(kyi)-tsamp])';

sprintf('Size of kxt %d,  nk %d',length(kxt),nk)

if flag
  nk = length(kxt)-2;
end


kxo = kxt(1:nk);
kyo = kyt(1:nk);


sprintf('Performing %d rotations',nint)
kx(:,1) = kxo;
ky(:,1) = kyo;
phi = 2*pi/nint;
for ii = 1:(nint-1)
     kx(:,ii+1) = kxo*cos(ii*phi) - kyo*sin(ii*phi);
     ky(:,ii+1) = kyo*cos(ii*phi) + kxo*sin(ii*phi);
end

kx = kx(:);
ky = ky(:);







