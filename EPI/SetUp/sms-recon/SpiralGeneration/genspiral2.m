function [Gx, Gy, kx, ky, sx, sy]=genspiral(D, N, Tmax, dts, nl); 

%   multi- shot spiral design 
%    uses Duyn's approximate slewrate limited design 
%    augmented with archimedian gmax limit 
%    inputs (args) 
%        D = FOV, cm 
%        N = matrix size 
%	 Tmax = longest acquisition allowed, s 
%	 dts = output sample spacing, s 
%        gtype = trajectory type 
%    inputs (CVs) 
%        nl = number of interleaves 
%        gamp = design grad max, G/cm 
%       gslew = design slew rate, mT/m/ms 
%	 nramp = number of rampdown points 
%    outputs 
%        Gx, Gy 
%        grev 
%    time is in sec 
% 
%		rev 0 12/26/98	original 
%		rev 1 4/15/99	little better calc of ts 
% 


%%%%%%%%%% Predefined variables

GRESMAX= 21000;
%nl=2;  %1   % Number of interleaves
gamp=2.2; %3.50; % 2.2 for both 1.5 T and 3 T data
gslew=120 % 200 % 180 for 3T data and 120 (150) for 1.5 T data
nramp=0;
% nramp=100;
MAX_PG_WAMP=32766;


gts = dts;
opfov = D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  gamma = 2*pi*4.257e3; 
gambar = gamma/(2*pi);

 
   gx=zeros(1,2*GRESMAX);
   gy=zeros(1,2*GRESMAX); 
 
 
  q = 5;
  S0 = gslew*100; 
  dt = dts*.5; 
 
%  slew-rate limited approximation  */ 
 
  Ts = .666667/nl*sqrt(power((pi*N), 3.0)/(gamma*D*S0)); 
  if(Ts > Tmax) display('slew limited readout too long'); end; 
   
  a2 = N*pi/(nl*power(Ts, .666667)); 
  a1 = 1.5*S0/a2; 
  beta = S0*gamma*D/nl; 
  Gmax = a1*power(Ts, .333333); 
  gmax = 0; 
  i = 1; 
  
  Ts
  
  for t = 0:dt:Ts, 
    x = power(t, 1.333333); 
    theta = .5*beta*t*t/(q + .5*beta/a2*x); 
    y = q+.5*beta/a2*x; 
    dthdt = beta*t*(q+.166667*beta/a2*x)/(y*y); 
    c = cos(theta); 
    s = sin(theta); 
    gx(i) = nl/(D*gamma)*dthdt*(c - theta*s); 
    gy(i) = nl/(D*gamma)*dthdt*(s + theta*c); 
    gabs = sqrt((gx(i)).^2 + (gy(i)).^2); 
    if(gabs>=gamp)   
       if(gmax==0)
          gmax = sqrt((gamp/theta).^2 + gamp.^2); 
       end;
       if(gabs>gmax)
         break;
       end;
    end;
    ts = t; 
    thetas = theta; 
    i=i+1; 
   end;
   
   
% gmax limited approximation  
 
  if(Gmax > gamp)   
    T = ((pi*N/nl)*(pi*N/nl) - thetas*thetas)/(2*gamma*gamp*D/nl) + ts; 
    if(T > Tmax) display('gmax limited readout too long'); end; 
     
    for t=(ts+dt):dt:T,
      theta = sqrt(thetas*thetas + 2*gamma*gamp*D*(t-ts)/nl); 
      c = cos(theta); 
      s = sin(theta); 
      gx(i) = gamp*(c/theta - s); 
      gy(i) = gamp*(s/theta + c); 
      i=i+1;
    end;
  end;
 
  %  decimate by 2 to get back to 4us sampling 
  
  n = 1; 
  totx=0;
  toty=0;
  
  for j=1:2:(i),   
    x = max(-gamp, min(gx(j), gamp)); 
    gxi = x*MAX_PG_WAMP/gamp; 
    Gx(n) = 2*(gxi/2); 
    y = max(-gamp, min(gy(j),gamp)); 
    gyi = y*MAX_PG_WAMP/gamp; 
    Gy(n) = 2*(gyi/2);
    
    totx=totx + Gx(n)*(gamp*gts*opfov*gambar/MAX_PG_WAMP);
    kx(n)= totx;
    toty=toty + Gy(n)*(gamp*gts*opfov*gambar/MAX_PG_WAMP);
    ky(n)= toty;
    
    n=n+1;
    
  end; 
  
  n
  
  
  
  bx = Gx(n-1); 
  by = Gy(n-1);
  
  for j=1:nramp,   
    c = 1 - (j-1)/nramp; 
    Gx(n) = 2*((bx*c)/2); 			% No int, problem?
    Gy(n) = 2*((by*c)/2);
     
    totx=totx + Gx(n)*(gamp*gts*opfov*gambar/MAX_PG_WAMP);
    kx(n)= totx;
    toty=toty + Gy(n)*(gamp*gts*opfov*gambar/MAX_PG_WAMP);
    ky(n)= toty;

    n=n+1;
  end; 
 
  gres = n - 1
 
  %kx=kx/(gamma/2*pi);
  %ky=ky/(gamma/2*pi);
  
  sx=diff(Gx)/dts;
  sx=sx(1:(length(sx)-1));
  
  sy=diff(Gy)/dts;
  sy=sy(1:(length(sy)-1));
 
