function ob = fast_mr_wt_quad(kx,ky,fov,N,K,J,tt,hanning_flag,we,r2,mask,flag_swt,L,int_opt,we_histo,epiflag,nx,ny,nz,zop,gz)
%function ob = fast_mr_wt(kx,ky,fov,N,K,J,tt,we,flag_swt,L,int_opt,we_histo,epiflag,nx,ny,nz,zop,gz)
%         usage such as A = fast_mr(kx,ky,FOV,N,2*N,5,tt,we(:),1,5,1);
%
%  or   ob = fast_mr(A_old,str_to_update,value)
%       usage such as A_new = fast_mr(A_old,'we',we_new(:))
%	Construct MRI object, which can do Ax and A'y operations
% typical input params:
%   kx goes form -N/2 to N/2:  unitless version
%   ky goes from -N/2 to N/2:  unitless version
%   fov is the fov (cm)
%    N = 64,128, MTX size
%    K oversampling   = 2*N, 3*N, 1.5*N
%    J = number of neighbors, typ. 5,6,7,8,...
%    flag_swt = 1 for sinc-weighting for rect basis functions
%    L number of time segments.   L=5,6,7,...
%    we_histo is a histogram of field map (we) for int_opt=2
%            has two columns:
%                    column 1:  bin centers for histogram
%                    column 2:  histogram values at those bin centers
%    int_opt = 1 for using input field map to find min-max int.
%    int_opt = 2 for histogram and include we_histo for histogram of field map
% Note:  if we = zeros, then L=0 is used.
%  Brad Sutton   Jeff Fessler
%  University of Michigan
%  October 28,2002
%  nx is the number oversampled in x
%  ny is the number oversampled in y
%  nz is the number of slices in z
%  zop determines whehter or not to inpterpolate input 3-D field map, or
%      whether a gz gradient is provided
%  zop = 1 mean use gradient
% This is an attempt to make MRI object valid for sensitivity
%   encoded runs also

%	default object
ob.G = 0;
ob.we = 0;
ob.tt = 0;
ob.int = 0;
ob.flgswt = 0;
ob.swt = 0;
ob.N = 0;
ob.is.empty	= logical(1);
ob.is.transpose = logical(0);
ob.r2=0;

if ~exist('zop','var')
  zop =0;
end


if nargin == 0
	ob = class(ob, 'fast_mr_wt_quad');
	return
end

if isa(kx, 'fast_mr_wt_quad')
	ob = kx;
        if nargin == 1
           return
        elseif (isstr(ky))
           eval(sprintf('ob.%s = fov;',ky))
        end
	return
end

if ((nargin < 18)| (nargin > 20))
	help fast_mr_wt
	error nargin
end



%	fill object
N1 = ny*N;
N2 = nx*N;
Jm1 = J;
Jm2 = J;
K1 = ny*K;
K2 = nx*K;
% ob.st = nufft2_init([2*pi/(N1)*ky,2*pi/(N2)*kx],N1,N2,Jm1,Jm2,K1,K2,[N/2-1/2,N/2-1/2]*(1),[],'kb');
%ob.st = nufft_init([2*pi/(N1)*ky,2*pi/(N2)*kx],[N1,N2],[Jm1,Jm2],[K1,K2],[N1/2-1/2,N2/2-1/2]*(1));
args = {[2*pi/N1*ky,2*pi/N2*kx],[N1,N2],[Jm1,Jm2],[K1,K2],[N1/2-1/2,N1/2-1/2]*(1)};
%args = {[2*pi/N1*ky,2*pi/N2*kx],[N1,N2],[Jm1,Jm2],[K1,K2],[N/2,N/2]*(1)};

if ~epiflag
  ob.G = Gnufft(args);
    else
  ob.G = Gfft(N1/nx,N1);
    if ~(K == N)
        sprintf('Oversmapling not used with epiflag')
    end
end


if (ndims(we) == 2)
we = reshape(we,N,N); % issue with jointest and four_echo_mri, field map is in column
r2 = reshape(r2,N,N);
end


H1 = sparse(kron(eye(N),1/(ny)*ones(ny,1)));
H2 = sparse(kron(eye(N),1/(nx)*ones(nx,1)));
HH = sparse(kron(H2,H1));

%Resize Field map
we_out = zeros(ny*N,nx*N,nz);
r2_out = zeros(ny*N,nx*N,nz);


% if (ndims(we) == 2)
%     nzwe = 1;
%     szx = size(we,2);
%     szy = size(we,1);
% else
    nzwe = size(we,3);
    szx = size(we,2);
    szy = size(we,1);
% end
wem = zeros(szy,szx,nz);

if zop == 1

    if ~exist('gz','var')
          sprintf('No z gradient detected')
    else
      gz = reshape(gz,szy,szx);
    end

    if ~(size(we,3)==1)
      sprintf('using zop = 1 requires fm to be single slice')
      keyboard
    end

    vout = [-nz/2:nz/2-1]./nz + 1/(2*nz);
      for ii = 1:szy
          for jj = 1:szx
             wem(ii,jj,:) = gz(ii,jj)*vout+we(ii,jj);
          end
      end
    we = wem;

else

    if ~(nzwe == nz)
      vin = [-nzwe/2:nzwe/2-1]./nzwe+1/(2*nzwe);
      vout = [-nz/2:nz/2-1]./nz + 1/(2*nz);

      for ii = 1:szy
          for jj = 1:szx
              P = polyfit(vin',squeeze(we(ii,jj,:)),1);
              wem(ii,jj,:) = P(1)*vout'+P(2);
          end
      end

    else
    wem = we;

    end
    we = wem;

end


for ii = 1:nz

   if ~(ndims(we) == 2)
       wec = squeeze(we(:,:,ii));
   else
       wec = we;
   end

   [nwey,nwex] = size(wec);

   if ~((nwey == ny*N)&&(nwex == nx*N))
       [wxl1,wyl1]= meshgrid([-(nwex)/2:(nwex)/2-1]./(nwex),[-(nwey)/2:(nwey)/2-1]./(nwey));
       [wxl2,wyl2]= meshgrid([-(nx*N)/2:(nx*N)/2-1]./(nx*N),[-(ny*N)/2:(ny*N)/2-1]./(ny*N));
       wen = interp2(wxl1+1/(2*nwex),wyl1+1/(2*nwey),wec,wxl2+(1/(2*nx*N)),wyl2+(1/(2*ny*N)),'cubic');
    r2n= interp2(wxl1+1/(2*nwex),wyl1+1/(2*nwey),r2,wxl2+(1/(2*nx*N)),wyl2+(1/(2*ny*N)),'cubic');
%       r2n = reshape(HH*r2(:),nx*N,ny*N);
   else
       wen = wec;
       r2n=r2;
   end

   % Smoothing the interpolated field map

%     beta = 1e3;
%     niter = 40;
%     init_img = zeros(nx*ny*N*N,1);
%     C = C2sparse('leak', ones(nx*N,ny*N), 8, 0);
%     C = sqrt(2^(-6)) * C;
%
%     Cn = Roughness_Penalty(ones(nx*N,ny*N))*beta;
%     I = ones(nx*N,ny*N);
%     A = spdiag(I(:));
%
%     se1 = strel('disk',2);
%      mask = reshape(HH*mask(:),nx*N,ny*N);
%     mask1 = imerode(mask,se1);
%     W = spdiag(mask1); % mask is the right size
%
%     tmp = qpwls_pcg(init_img, A, W,col(wen), 0, C, 1, niter);
%     tmp = reshape(tmp,nx*N,ny*N,niter);
%    r2n = r2; %interpolation of r2 done in mtimes using HH matrice
%    we_out(:,:,ii) = tmp(:,:,end);
   we_out(:,:,ii) = wen;
   r2_out(:,:,ii) =r2n;

end



llwe = find(isnan(we_out));
we_out(llwe) = 0;

llr2 = find(isnan(r2_out));
r2_out(llr2) = 0;


%  ob.int = timr2_seg_int_mex_mex(tt,L,we(:),r2(:),3, 1000);
if hanning_flag
    ob.int = int_tim_seghanning(tt,L);
else
    ob.int = timr2_seg_int(tt,L,we(:),r2(:),3, 1000);

end
% AA = int_tim_seg(tt,L,we(:),1, []);
ob.N = N;

if flag_swt
   ob.flgswt = 1;
   ob.swt = sinc(kx/N2).*sinc(ky/N1);
   %ob.swt = sinc(kx/N).*sinc(ky/N);
else
   ob.flgswt = 0;
end


ob.we = we_out;
ob.tt = tt;
ob.is.empty	= logical(0);

%set up HH

ob.HH = HH;

ob.r2=r2_out;

ob = class(ob, 'fast_mr_wt_quad');
