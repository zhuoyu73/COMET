% This file will be used to build a database for the generic interpolator

dest = '/net/waters/home/bpsutton/Matlab/@fast_mr/Int';

FOVlist = [20;22;24];
nllist = [1;2;4;8;16];

NN = 512;
N = 64;
gamp = 2.2;
gslew = 180;
rotamount = 0;
rev_flag = 0;
sh = 0;
sh_pos = 1;
TE = 25e-03;

for indfov = 1:length(FOVlist)
    
    for indnl = 1:length(nllist)
FOV = FOVlist(indfov);

nl = nllist(indnl);
L = 8;
  % %K-SPACE TRAJECTORY
%**********************************************************
% Load in k-space trajectory

if 1
  Tmax = 5*16384*1e-6;   %5*16384*1e-6;
  dts = 5e-6;  %info3(14)*1e-06;    %5e-6
  [kx,ky] = genkspace(FOV,N,0,nl,gamp,gslew,dts,rotamount,rev_flag);
  ndat = length(kx)/nl;    %num data pts for one spiral
  datseg = ndat;          %set num data pts to one spiral
  
  %shifting k-space data
  for ll = 1:nl
    if sh_pos
      kx_set = sprintf('kx%d = kx(((ll-1)*ndat+1):(ll*ndat-sh));',ll);
      eval(kx_set);
      ky_set = sprintf('ky%d = ky(((ll-1)*ndat+1):(ll*ndat-sh));',ll);
      eval(ky_set);
    else
      kx_set = sprintf('kx%d = kx(((ll-1)*ndat+1+sh):(ll*ndat));',ll);
      eval(kx_set);
      ky_set = sprintf('ky%d = ky(((ll-1)*ndat+1+sh):(ll*ndat));',ll);
      eval(ky_set);  
    end
  end
  %Loading kx1,... back into kx
  kx_setup = 'kx = [kx1';
  ky_setup = 'ky = [ky1';
  for mm = 1:nl-1
    kx_setup = sprintf('%s; kx%d',kx_setup,mm+1);
    ky_setup = sprintf('%s; ky%d',ky_setup,mm+1);
  end
  kx_setup = sprintf('%s];',kx_setup);
  ky_setup = sprintf('%s];',ky_setup);
  eval(kx_setup);
  eval(ky_setup);
end

  if rev_flag
    ww = flipud(weight_vor(flipud(kx),flipud(ky),nl));
  else
    ww = weight_vor(kx,ky,nl);
  end
  
  


%TIMING
%********************************************************
% Set up timing of slice acquisition
if 1
  tt = [(0:(ndat-1))*dts+TE]'; %Timing of non-delayed acquisition
  %Shifting timing vector
  if sh_pos
    tt = tt((1+sh):ndat);
  else
    tt = tt(1:(ndat-sh));  
  end
  %Put back into timing vector for all interleaves, tt_ext
  tt_setup = 'tt_ext = [tt';
  for mm = 1:nl-1
    tt_setup = sprintf('%s; tt',tt_setup);
  end
  tt_setup = sprintf('%s];',tt_setup);
  eval(tt_setup);
end
 

bin_cens = linspace(-30*2*pi,160*2*pi,NN);
  bin_vals = (N*N)/NN*ones(NN,1);
   we_histo = [bin_cens', bin_vals];

   we = zeros(N*N,1);
   

AA = int_tim_seg(tt_ext,L,we,2, we_histo);
eval(sprintf('intl%dm%dnl%d = AA;',L,ndat,nl))
save(sprintf('%s/intl%dm%dnl%d',dest,L,ndat,nl),sprintf('intl%dm%dnl%d',L,ndat,nl))
eval(sprintf('fid = fopen(''%s/intl%dm%dnl%d'',''w'');',dest,L,ndat,nl))
fwrite(fid,[real(AA(:)) imag(AA(:))],'float')
fclose(fid)

end

end

