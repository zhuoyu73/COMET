function [kx,ky,gx,gy] = genkspace(FOV,N,ld,nint,gamp,gslew,tsamp,rotamount,rev_flag,gts,flag_vd,int_rotation,alpha_vd);
%[kx,ky,gxt,gyt] = genkspace(FOV,N,ld,nint,gamp,gslew,tsamp,rotamount,rev_flag,gts,flag_vd,int_rotation,alpha_vd);
%function [kxo,kyo] = genkspace(FOV,N,Tmax,ld,nint,gamp,gslew,gts);
% This function will generate the proper length of k-space
% trajectory.  It linearly interpolates the output of
% genspiral to the correct length and takes care of the 
% rotations for the interleaves.
% ld is the length of the data
% nint is the number of interleaves
% flag_vd is a switch for vairable density
% int_rotation is the rotation over interleave (given as number of interleaves)
%Brad Sutton, University of Michigan
flag = 0;   %auto determine number of k-space points
            % just input ld = 0.3

if ~(exist('rotamount','var'))
  rotamount = 0;
end
if ~(exist('rev_flag','var'))
  rev_flag = 0;
end
if ~(exist('gts','var'))
  gts = 4e-6;
end
if ~(exist('flag_vd','var'))
  flag_vd = 0;
end
if ~(exist('int_rotation','var'))
  int_rotation = 1;
end
if ~(exist('alpha_vd','var'))
  alpha_vd = 4;
end
if (int_rotation < 0)
  flag_sign = -1;
  int_rotation = -int_rotation;
else
  flag_sign = 1;
end

nk = ld/nint;
if round(nk) ~= nk
  sprintf('Input should have num data pts/number of interleaves must be int')
end

if (nk == 0)
     flag = 1;
end

%dts = 4e-6;    %5e-6
if (flag_vd == 1)
  [Gx,Gy,kxi,kyi,sx,sy] = genspivd_Kim(FOV,N,nint,gamp,gslew,gts,alpha_vd);
else
  [Gx,Gy,kxi,kyi,sx,sy] = genspi(FOV,N,nint,gamp,gslew,gts);
end

kxt=interp1([0:gts:gts*length(kxi)-gts],kxi,[0:tsamp:gts*length(kxi)-tsamp])';
kyt=interp1([0:gts:gts*length(kyi)-gts],kyi,[0:tsamp:gts*length(kyi)-tsamp])';

if nargout>2
  gxt=Gx;
  gyt=Gy;
end


if flag
     nk = length(kxt)-2;   
end

kx = zeros(nk,nint);
ky = zeros(nk,nint);
kxo = zeros(nk,1);
kyo = zeros(nk,1);


kxo = kxt(1:nk);
kyo = kyt(1:nk);


if nargout>2
  ng = length(gxt);
  gx = zeros(ng,nint);
  gy = zeros(ng,nint);
  gxo = gxt(:);
  gyo = gyt(:);

end

%if length(kxi)==nk;
%  kxo = kxi.';
%  kyo = kyi.';
%else

%  kxo(1) = kxi(1);
%  kyo(1) = kyi(1);

%  nki = length(kxi);
%  sprintf('Interpolating %d kspace pts to %d',nki,nk)

%  sf = (nki-1)/(nk-1);  %scaling factor
%  for ii = 1:(nk-2);
%    ind = ii*sf;
%    kxo(ii+1) = kxi(floor(ind)+1)*(1+floor(ind)-ind)+kxi(ceil(ind)+1)*(ind-floor(ind));
%    kyo(ii+1) = kyi(floor(ind)+1)*(1+floor(ind)-ind)+kyi(ceil(ind)+1)*(ind-floor(ind));
%  end
%  kxo(end) = kxi(end);
%  kyo(end) = kyi(end);
%end

%rotate matrix for proper orientation
phir = -rotamount*pi/2;
kxop = kxo*cos(phir) - kyo*sin(phir);
kyop = kyo*cos(phir) + kxo*sin(phir);


if rev_flag
  kxop = -flipud(kxop);
  kyop = -flipud(kyop);
end


sprintf('Performing %d rotations',nint)
%kx(:,1) = kxop;
%ky(:,1) = kyop;
phi = 2*pi/nint;

for ii = 0:(nint-1)
     ang_rot = phi*flag_sign*(ii*int_rotation - (nint-1)*floor(ii*int_rotation/nint));
     %ang_rot = flag_sign*(mod(int_rotation*(ii),nint)+floor(int_rotation*(ii)/nint))*phi;
     kx(:,ii+1) = kxop*cos(ang_rot) + kyop*sin(ang_rot);
     ky(:,ii+1) = kyop*cos(ang_rot) - kxop*sin(ang_rot);
end

kx = kx(:);
ky = ky(:);


if nargout>2

gxop = gxo*cos(phir) - gyo*sin(phir);
gyop = gyo*cos(phir) + gxo*sin(phir);
  %gx(:,1) = gxop;
  %gy(:,1) = gyop;
  for ii = 0:(nint-1)
     ang_rot = phi*flag_sign*(ii*int_rotation - (nint-1)*floor(ii*int_rotation/nint));
     %ang_rot = flag_sign*(mod(int_rotation*(ii),nint)+floor(int_rotation*(ii)/nint))*phi;
     gx(:,ii+1) = gxop*cos(ang_rot) + gyop*sin(ang_rot);
     gy(:,ii+1) = gyop*cos(ang_rot) - gxop*sin(ang_rot);
  end

  gx = gx(:);
  gy = gy(:);
end


%keyboard

