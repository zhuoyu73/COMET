function [kx,ky] = genspivd(opfov,opxres, nl, gamp, gslew, densamp)
%function [kx,ky] = genspivd(opfov,opxres, nl, gamp, gslew, densamp)
% From Doug Noll
%  University of Michigan
%  Borrowed on 11/29/01
%  see ~dnoll/epic/lx/spirals/dogradvd.m

if 0
   gslew = 180;		% mt/m/ms
   gamp = 2.2;		% fullscale g/cm

   opxres = 96;
   nl = 2;			% number of interleaves
   opfov = 24;		% cm
   densamp = 300;		% duration of higher density sampling
end

fsgcm = gamp;           % fullscale g/cm
risetime = gamp/gslew*10000;    % us
targetk = opxres/2;

ts = 4e-6;		% sampling time
gts = 4e-6;		% gradient sampling time
dentrans = densamp/2;	% duration of transition from higher to lower
A = 32766;		% output scaling of waveform (fullscale)

MAXDECRATIO = 32;
GAM = 4257.0;
S = (gts/1e-6)*A/risetime;
dr = ts/gts;
OMF = 2.0*pi * opfov/(1/(GAM*fsgcm*gts));
OM = 2.0*pi/nl * opfov/(1/(GAM*fsgcm*gts));
distance = 1.0 / (opfov*GAM*fsgcm*gts/A);
i = sqrt(-1);

ac = A;
loop = 1;
decratio = 1;
absk = 0;
S0 = gslew*100; 
ggx = [];
ggy = [];
dens = [];

while loop > 0,		% start over
	loop = 0;
	om = OM/decratio; 
	omf = OMF/decratio; 
	s = S/decratio;
	g0 = 0; 
	gx = g0; 
	gy = 0; 
	absg = abs(g0);
	oldkx = 0; 
	oldky = 0; 
	tkx = gx; 
	tky = gy; 
	kxt = tkx; 
	kyt = tky;
	thetan_1 = 0; 
	taun = 0; 
	n = 0;
	den1 = 0;
	while absk < targetk
	    realn = n/decratio;
	    if rem(realn,200) == 1
	      realn
	    end
	    taun_1 = taun; 
	    taun = abs(tkx + i*tky)/A; 
	    tauhat = taun;
	    if realn > densamp
	      if den1 == 0 
		kdenrad = abs(tkx + i*tky)/distance/decratio
		den1 = 1;
	      end
	      if realn > densamp+dentrans
		scoffset = scthat;
		denoffset = taun_1;
	        scthat = scoffset + om*(tauhat - denoffset);
		fractrans = 1;
	      else
		scoffset = scthat; 
		denoffset = taun_1;
		fractrans = (realn - densamp)/dentrans;
		fractrans = 1 - ( (fractrans-1)*(fractrans-1));
	        scthat = scoffset + (omf + (om-omf)*fractrans)*(tauhat - denoffset);
	      end
	    else
 	      fractrans = 0;
	      scthat = omf*tauhat;
	    end
	    theta = atan2(scthat,1.0)+scthat;
	    if absg < ac
		deltheta = theta-thetan_1;
		B = 1.0/(1.0+tan(deltheta)*tan(deltheta));
		gtilde = absg;
		t1 = s*s;
		t2 = gtilde*gtilde*(1-B);
		if t2 > t1
		    decratio = decratio * 2.0
		    if decratio > MAXDECRATIO
			printf('k-space calculation failed.\n');
			return;
		    end
		    loop = 1;
		    break;
		end
		t3 = sqrt(t1-t2);
		absg = sqrt(B)*gtilde+t3;
		if (absg > ac)
		    absg = ac;
	        end
	    end
	    tgx = absg*cos(theta);
	    tgy = absg*sin(theta);
	    tkx = tkx + tgx;
	    tky = tky + tgy;
	    thetan_1=theta;
	    if rem(n,decratio) == 0
		m = round(n/decratio);
		gx = round(((tkx-oldkx))/decratio);
		gx = gx - rem(gx,2);
		gy = round(((tky-oldky))/decratio);
		gy = gy - rem(gy,2);
		ggx(m+1) = gx;
		ggy(m+1) = gy;
		kxt = kxt + gx;
		kyt = kyt + gy;
		oldkx = tkx;
		oldky = tky;
	    	if rem(m,dr) == 0
		    m  = m / dr;
		    absk = abs(kxt + i*kyt)/distance;
		    dens(m+1) = omf/(omf + (om-omf)*fractrans);
		    if absk > targetk
			break;
		    end
		    kx(m+1) = kxt/distance;
		    ky(m+1) = kyt/distance;
		end
	    end
	    n = n+1;
	end
end
npts = m+1
arraysize = (2*nl*om*taun/(2*pi))
res = opfov/arraysize*10
gvec = (ggx + i.*ggy)./A.*fsgcm;
svec = diff(gvec)/gts/1000;
plot((1:npts-1).*gts.*1000,abs(svec),(1:npts).*gts.*1000,abs(gvec))

return

% calculate true k locations - 3 term integrations - simpson's rule
%   (simpson's rule 1/6 4/6 1/6 weighting)
dt = gts*1000;
delk = 1/4.258/opfov; 	% (g ms)/cm
k(1) = 0;
k(2) = (gvec(2) + 4*gvec(1) )*dt/6;
for i=1:npts-2,
  k(i+2) = k(i+1) + (gvec(i+2) + 4*gvec(i+1) + gvec(i))*dt/6;
end % for i=1:npts-1
% create .k file
outk = k(1:npts) ./delk;

return

filename = sprintf('s%df%d.k',nl,opfov);
fid = fopen(filename,'w');
for lp=0:(nl-1),
  fwrite(fid,npts,'int');
  fwrite(fid,real(outk.*exp(sqrt(-1).*(2.*pi.*lp./nl+ pi/2))),'float');
  fwrite(fid,imag(outk.*exp(sqrt(-1).*(2.*pi.*lp./nl+ pi/2))),'float');
end % for i=1:nl,
fclose(fid);

