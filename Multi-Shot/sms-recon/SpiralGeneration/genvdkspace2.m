function [kx,ky] = genvdkspace2(direc,kfile,ld,nint);
%function [kx,ky] = genvdkspace2(filename,ld,nint);
% This function will generate the proper length of k-space
% trajectory.  It linearly interpolates the output of
% genspiral to the correct length and takes care of the 
% rotations for the interleaves.
% ld is the length of the data
% nint is the number of interleaves
% filename is the location of the filename with traj in it.

nk = ld/nint;
if round(nk) ~= nk
sprintf('Input should have num data pts/number of interleaves must be int')
end

kxo = zeros(nk,1);
kyo = zeros(nk,1);
kx = zeros(nk,nint);
ky = zeros(nk,nint);

dts = 4e-6;    %5e-6
eval(sprintf('load %s/%s',direc,kfile));
eval(sprintf('kxi = %s(:,1);',kfile));
eval(sprintf('kyi = %s(:,2);',kfile));

kxi = [0;kxi];
kyi = [0;kyi];

if length(kxi)==nk;
  kxo = kxi;
  kyo = kyi;
else

  kxo(1) = kxi(1);
  kyo(1) = kyi(1);

  nki = length(kxi);
  sprintf('Interpolating %d kspace pts to %d',nki,nk)

  sf = (nki-1)/(nk-1);  %scaling factor
  for ii = 1:(nk-2);
    ind = ii*sf;
    kxo(ii+1) = kxi(floor(ind)+1)*(1+floor(ind)-ind)+kxi(ceil(ind)+1)*(ind-floor(ind));
    kyo(ii+1) = kyi(floor(ind)+1)*(1+floor(ind)-ind)+kyi(ceil(ind)+1)*(ind-floor(ind));
  end
  kxo(end) = kxi(end);
  kyo(end) = kyi(end);
end


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







