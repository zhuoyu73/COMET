function [Gx, Gy, kx, ky, sx, sy]=genspivd(D, N, Tmax, dts, nl, gamp, gslew, dsamp, sfact); 
int getrttrajvd (opxres, nl, dsamp, sfact, dts, gts, gamp, opfov, gslew, Tmax, gtype, kx, ky) 

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

%NOTE TO SELF, Tmax = DTSMAX*5e-06?


%%%%%%%%%% Predefined variables

GRESMAX= 21000;
if ~exist('nl','var')
    nl=1   % Number of interleaves
end
if ~exist('gamp','var')
      gamp=2.2; %3.50; % 2.2 for both 1.5 T and 3 T data
end
if ~exist('gslew','var')
      gslew=180 % 200 % 180 for 3T data and 120 (150) for 1.5 T data
end
nramp=0;
% nramp=100;
MAX_PG_WAMP=32766;
MAXDECRATIO = 32;

gts = dts
opfov = D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  gamma = 2*pi*4.257e3; 
gambar = gamma/(2*pi);

%JUST INCLUDED ~dnoll/recon/getrttrajvd.c and modified to run in Matlab

MAXDECRATIO  = 32 % maximum allowed decimation of input ts 
MAX_PG_WAMP = 32766
GAM = 4257.0       
GRESMAX = 21000
DTSMAX = 16384


    r1 = zeros(nl,1);
    r2 = zeros(nl,1);
    for ii = 1:nl
%	r1(ii) = cos(2*M_PI*((double)i/nl+0.25));
%	r2(ii) = sin(2*M_PI*((double)i/nl+0.25));
	r1(ii) = cos(2*pi*(ii/nl+0.25));
	r2(ii) = sin(2*pi*(ii/nl+0.25));
    end

kx = zeros(nl,DTSMAX);
ky = zeros(nl,DTSMAX);


    A = MAX_PG_WAMP;
    risetime = gamp/gslew*10000;
    S = (dts/1e-6)*A/risetime;

    targetk = opxres/2;

    A = MAX_PG_WAMP;
    S = (gts/1e-6)*A/risetime;
%    OM = (2*M_PI)/(sfact*nl) * (opfov)/(1/(GAM*gamp*gts));
%    OM = (2*pi)/(sfact*nl) * (opfov)/(1/(GAM*gamp*gts));

    /* distance for one k-space unit */
    distance = 1.0 / (opfov*GAM*gamp*gts/A); 
    OMF = OM*sfact*nl;
    dentrans = dsamp/2;

    /* printf("OM = %f, S = %f, A = %f\n", OM, S, A); */
    ac = A; loop = 1; decratio = 1;
    while (loop)
        absk = 0.;
	loop = 0;
	om = OM/decratio;
	s = S/decratio;
	omf = OMF/decratio;
	den1 = 0;
	g0 = 0;
	absg = hypot(g0,0);
	oldkx = 0;
	oldky = 0;
	tkx = g0;
	tky = 0;
	kxt = 0;
	kyt = 0;
	thetan_1 = 0;
	taun = 0;
	n = 0;
	while (absk < targetk)
	    taun_1 = taun;
	    taun = hypot(tkx,tky)/A;
	    tauhat = taun;
	    realn = ((float) n)/((float) decratio);
	    if (realn > (float) dsamp) {
	      if (den1 == 0)  {
		ksv = taun;
		den1 = 1;
		/* printf("kdenrad = %f\n", sfact*nl*om*ksv/(2*M_PI)); */
	      }
	      if (realn > (float) (dsamp+dentrans)) {
		scoffset = scthat;
		denoffset = taun_1;
	        scthat = scoffset + om*(tauhat - denoffset);
		fractrans = 1;
	      }
	      else {
		scoffset = scthat;
		denoffset = taun_1;
		fractrans = (realn - dsamp)/((float) dentrans); 
		fractrans = 1 - ( (fractrans-1)*(fractrans-1));
	        scthat = scoffset + (omf + (om-omf)*fractrans)*(tauhat - denoffset);
	      }
	    }
	    else {
	      scthat = omf*tauhat;
	      fractrans = 0;
	    }
	    theta = atan2(scthat,1.0)+scthat;
	    if (absg < ac)
	    {
		deltheta = theta-thetan_1;
		B = 1.0/(1.0+tan(deltheta)*tan(deltheta));
		gtilde = absg;
		t1 = s*s;
		t2 = gtilde*gtilde*(1-B);
		if (t2 > t1)
		{
		    decratio *= 2.0;
			/* printf("decratio = %f\n",decratio); */
		    if (decratio > MAXDECRATIO)
		    {
			printf("gettrajrtvd failed\n");
			      exit(0);
		    }
		    loop = 1;
		    break;
		}
		t3 = sqrt(t1-t2);
		absg = sqrt(B)*gtilde+t3;
		if (absg > ac)
		    absg = ac;
	    }
	    tgx = absg*cos(theta);
	    tgy = absg*sin(theta);
	    tkx += tgx;
	    tky += tgy;
	    thetan_1=theta;
	    if (!(n % ((int) rint(decratio))))
	    {
		m = n/((int) rint(decratio));
		gx = ((int) rint((tkx-oldkx))/decratio);
		gy = ((int) rint((tky-oldky))/decratio);
		kyt += gy;
		absk = hypot(kxt,kyt)/distance;
		oldkx = tkx;
		oldky = tky;
		Gx1[m] = gx;
		Gy1[m] = gy;
		/* if (!(m%200))
		  printf("m = %d, absk = %f, gx,gy = %d, %d\n",m,absk,gx,gy); */
	    }
	    n++;
	end %end while
    end  %end while
  gres = m + 1; 
  Kx1[0] = 0.;
  Ky1[0] = 0.;
  for (i = 0, kxt = 0, kyt= 0; i < gres; i++) {
    gx = Gx1[i];
    gy = Gy1[i];
    gx &= 0xfffe;
    gy &= 0xfffe; 
    kxt += gx;
    kyt += gy;
    Kx1[i+1] = kxt/distance;
    Ky1[i+1] = kyt/distance;
    absk = hypot(kxt,kyt)/distance;
  }
    res = ceil(2*sfact*nl*om*taun/(2*M_PI));
    kdenrad = ceil(2*sfact*nl*om*ksv/(2*M_PI));
    if (out_flg)
  printf("rad = %f, high dens region = %f, npts = %d\n", absk, (2*sfact*nl*om*ksv/(2*M_PI)),gres); 

   
  for (i = 0; i < nl; i++)
  {
    kx[i][0] = 0.0;
    ky[i][0] = 0.0;
  }
  n=0;
  t = gres*gts;
  for (ts=dts; ts<=t; ts+=dts)  { 
    indl = floor(ts/gts);
    if (indl < gres) {
      kxt = Kx1[indl]*(1- ts/gts + indl) + Kx1[indl+1]*(ts/gts - indl);
      kyt = Ky1[indl]*(1- ts/gts + indl) + Ky1[indl+1]*(ts/gts - indl);
      for (i = 0; i < nl; i++)
      {
       kx[i][n+1] = (kxt*r1[i] - kyt*r2[i]);
       ky[i][n+1] = (kxt*r2[i] + kyt*r1[i]);
      }
      n++;
    }
  }
  return n;
}    

float *fspace2(size)
int size;
{
    float *buffer;
    if (!(buffer = (float *)calloc( size,sizeof (float) ) ))
    { fprintf(stderr,"calloc: space not assigned"); exit (0); }
    return(buffer);
}








XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXgenspiral code here 
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
%    totx=totx + Gx(n);
    kx(n)= totx;
    toty=toty + Gy(n)*(gamp*gts*opfov*gambar/MAX_PG_WAMP);
%    toty=toty + Gy(n);
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
%    totx=totx + Gx(n);
    kx(n)= totx;
    toty=toty + Gy(n)*(gamp*gts*opfov*gambar/MAX_PG_WAMP);
%    toty=toty + Gy(n);
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
 
