function [Gx, Gy, kx, ky, sx, sy]=genspi(D, N, nl, gamp, gslew, gts); 
%function [Gx, Gy, kx, ky, sx, sy]=genspi(D, N, nl, gamp, gslew, gts);
%   multi- shot spiral design 
%    uses Duyn's approximate slewrate limited design 
%    augmented with archimedian gmax limit 
%    inputs (args) 
%        D = FOV, cm 
%        N = matrix size
%	 %Tmax = longest acquisition allowed, s 
%	 gts = output sample spacing, s, gradient raster time 
%        gtype = trajectory type 
%    inputs (CVs) 
%        nl = number of interleaves 
%        gamp = design grad max, G/cm 
%        gslew = design slew rate, mT/m/ms 
%	 nramp = number of rampdown points 
%    outputs 
%        Gx, Gy 
%        grev 
%    time is in sec 
% 
%		rev 0 12/26/98	original 
%		rev 1 4/15/99	little better calc of ts 
% 
%               borrowed from Doug Noll, Univ. of Michigan
%               modified to take more input cv's 
% IMplementation described in Glover, Simple Analytic spiral k-space
% algorithm  MRM 42:412-415.
  
%%%%%%%%%% Predefined variables

GRESMAX= 32768;  % Value from GE, Siemens value??
if ~exist('nl','var')
    nl=1   % Number of interleaves
end
if ~exist('gamp','var')
      gamp=4; %2.2 for GE
      		% 4 for Allegra
end
if ~exist('gslew','var')
      gslew=400; % 180 for GE
      		% 400 for Allegra
end
if ~exist('gts','var')
  gts = 10e-06;   % 4e-06 for GE
  		 %10e-06 for Siemens
end
  nramp=0;
% nramp=100;
MAX_PG_WAMP=32766;

Tmax = GRESMAX*gts;

dts = gts;
opfov = D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma = 2*pi*4.257e3; 
gambar = gamma/(2*pi);

 
   gx=zeros(1,2*GRESMAX);
   gy=zeros(1,2*GRESMAX); 
 
 
  q = 5;
  S0 = gslew*100; 
  dt = dts*.5; 
 
%  slew-rate limited approximation  */ 
 
  Ts = .666667/nl*sqrt(((pi*N)^3)/(gamma*D*S0)); 
  if(Ts > Tmax) display('slew limited readout too long'); return; end; 
   
  a2 = N*pi/(nl*(Ts^(.666667))); 
  a1 = 1.5*S0/a2; 
  beta = S0*gamma*D/nl; 
  Gmax = a1*(Ts^.333333); 
  gmax = 0; 
   
  Ts
  
    t = [0:dt:Ts]; 
    x = t.^1.333333; 
    theta = (t.^2).*(.5*beta./(q + .5*beta./a2.*x)); 
    y = q+.5.*beta./a2.*x; 
    dthdt = t.*(beta.*(q+.166667*beta./a2.*x)./(y.*y)); 
    c = cos(theta); 
    s = sin(theta); 
    gx = (nl/(D*gamma)).*dthdt.*(c - theta.*s); 
    gy = (nl/(D*gamma)).*dthdt.*(s + theta.*c); 
    gabs = abs(gx+i.*gy); 
% cut short if over peak
  gmax = abs(gamp./(theta+eps) + i.*gamp); 
  l1 = length(t) - sum(gabs>gmax);
  ts = t(l1);
  thetas = theta(l1);

 
%/*  gmax limited approximation  */ 
 
l3 = 0;
T=ts;
if Gmax > gamp 
  T=((pi*N/nl)*(pi*N/nl) - thetas*thetas)/(2*gamma*gamp*D/nl)+ts;
  if T > Tmax 
      sprintf('gmax limited readout too long')
      return;
  end;
  t = [ts+dt:dt:T];
  theta = sqrt(thetas*thetas + (2*gamma*gamp*D).*(t-ts)./nl); 
  c = cos(theta); 
  s = sin(theta); 
  ind2 = l1+[1:length(t)];
  gx(ind2) = gamp.*(c./theta - s);
  gy(ind2) = gamp.*(s./theta + c);
  l3 = length(t);
end; 
 
l2 = l1 + l3;
Gx = gx(1:2:l2);   % or gx(1:2:l2)*MAX_PG_WAMP/gamp
Gy = gy(1:2:l2);   % or gy(1:2:l2)*MAX_PG_WAMP/gamp
g = Gx + i.*Gy;   %slew rate vector
s = diff(g)./(gts*1000);  % grad vector
Kx = cumsum([0 Gx])*gts*opfov*gambar;
Ky = cumsum([0 Gy])*gts*opfov*gambar;
k = Kx + i.*Ky;  %kspace vector
t = [0:gts:T]; % time vector
matrix = max(abs(k))*2
maxg = max(abs(g))
maxs = max(abs(s))
maxt = max(t).*1000

kx = real(k);
ky = imag(k);
sx = real(s);
sy = imag(s);
