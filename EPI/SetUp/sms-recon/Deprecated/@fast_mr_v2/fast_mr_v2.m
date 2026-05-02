function ob = fast_mr_v2(kx,ky,kz,fov,Nx, Ny, Nz,Kx,Ky,Kz, J,tt,we,flag_swt,L,int_opt,we_histo, gridflag,mask)
%function ob = fast_mr(kx,ky,fov,N,K,J,tt,we,flag_swt,L,int_opt,we_histo,epiflag)
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
%    int_opt = 2 for histogram and include we_histo for histogram of
%    field map 
%    gridflag = 0 (default), will performing grdding with 2 times OS
%    gridflag = 1 (EPI)  No need to oversample, just fft.
%    gridflag = 2 will performing grdding with 2 times OS, but use the
%               table look-up for memory efficiancy


% Note:  if we = zeros, then L=0 is used.
%
%  Brad Sutton   Jeff Fessler 
%  University of Michigan
%  June 17, 2002 
%
% History:  added epiflag 4/9/04
%           added mask    12/16/2013 JLH
%
% This is an attempt to make MRI object valid for sensitivity
%   encoded runs also
%
%..............................ANH's EDITS.................................
%
% Purpose: to speed up in the case of multi-shot, multi kz lines per
% volume
%    - nl: number of shot, to save computation time when performing
%          field map correction on a multishot acquisition
%    - nkz: number of kz encoding lines per volume
%    - tt: timing vector for ONE shot (in a multi-shot) acquisition
%
%
%..............................Joe's EDITS.................................
%   - Changed tt so that it now is a timing vector for each data point
%       This allows for arbitrary readouts, all readouts no longer need to
%       be the same
%   - Added the mask input, this allows for modeling the extent of the
%       object
%   - The EPI flag is now a more general grdding options flag
%   -tt can now be the timing vector for all the data or just for a single
%       shot (assuming all shots have the same tt)
warning('Deprecated: This object is deprecated and not recommended for new reconstruction routines. Use NUFFT or TimeSegmentation instead.');
%	default object
ob.G = 0;
ob.we = 0;
ob.tt = 0;
ob.int = 0;
ob.flgswt = 0;
ob.swt = 0;
ob.is.empty	= logical(1);
ob.is.transpose = logical(0);
ob.is_masked = 0;
ob.mask = [];
ob.ttsegments = 1;
%ob.version = 1.0;

if nargin == 0
	ob = class(ob, 'fast_mr_v2');
	return
end


if isa(kx, 'fast_mr_new')
	ob = kx;
        if nargin == 1
           return
        elseif (isstr(ky))
           eval(sprintf('ob.%s = fov;',ky)) 
        end
	return
end

if (nargin < 16)|(nargin>22)
	help fast_mr_v2
	error nargin
end

if ~exist('gridflag','var')
  gridflag = 0;
end

if ~exist('int_opt','var')
  int_opt = 1;
end


% if ~exist('mask','var')
%   mask = true(Nx,Ny,Nz);
% end

%	fill object
N1 = Nx;
N2 = Ny;
N3 = Nz;
Jm1 = J;
Jm2 = J;
if Nz == 1
    Jm3 = 1;
elseif Nz < 8
    Jm3 = 4;
else
    Jm3 = J;
end
K1 = Kx; 
K2 = Ky;
K3 = Kz;

%ob.st =
%nufft2_init([2*pi/N1*ky,2*pi/N2*kx],N1,N2,Jm1,Jm2,K1,K2,[N/2-1/2,N/2-1/2]*(1),[],'kb');
if(length(kz) > 1)
    warning('Performing 3D NUFFT...')
    if gridflag == 2
        args = {[2*pi/N1*ky,2*pi/N2*kx, 2*pi/N3*kz],[N1,N2, N3],[Jm1,Jm2,Jm3],[K1,K2,K3],[Nx/2-1/2,Ny/2-1/2, Nz/2-1/2], 'table', 2^10, 'minmax:kb'}; 
        warning('Using table lookup...');
    else
        args = {[2*pi/N1*ky,2*pi/N2*kx, 2*pi/N3*kz],[N1,N2, N3],[Jm1,Jm2,Jm3],[K1,K2,K3],[Nx/2-1/2,Ny/2-1/2, Nz/2-1/2], 'minmax:kb', 2^10}; 
    end
%     args = [ {[2*pi/N1*ky,2*pi/N2*kx, 2*pi/N3*kz],[N1,N2, N3],[Jm1,Jm2,Jm3],[K1,K2,K3]}, args]; 
%     args = {[2*pi/N1*ky,2*pi/N2*kx, 2*pi/N3*kz],[N1,N2, N3],[Jm1,Jm2,Jm3],[K1,K2,K3],[Nx/2,Ny/2, Nz/2], 'table', 2^10, 'minmax:kb'}; 
%     args = {[2*pi/N1*ky,2*pi/N2*kx, 2*pi/N3*kz],[N1,N2, N3],[Jm1,Jm2,Jm3],[K1,K2,K3],[Nx/2,Ny/2, Nz/2], 'minmax:kb', 2^10}; 
%     args = {[2*pi/N1*ky,2*pi/N2*kx, 2*pi/N3*kz],[N1,N2, N3],[Jm1,Jm2,Jm3],[K1,K2,K3],[Nx/2-1/2,Ny/2-1/2, Nz/2-1/2], 'linear', 2^10}; 
else
    if gridflag == 2
        args = {[2*pi/N1*ky,2*pi/N2*kx],[N1,N2],[Jm1,Jm2],[K1,K2],[N1/2-1/2,N2/2-1/2]*(1), 'table', 2^10, 'minmax:kb'};
        warning('Using table lookup...');
    else
        args = {[2*pi/N1*ky,2*pi/N2*kx],[N1,N2],[Jm1,Jm2],[K1,K2],[N1/2-1/2,N2/2-1/2]*(1), 'minmax:kb', 2^10};
    end
    
end

if (gridflag == 0)||(gridflag == 2)
    if ~exist('mask','var')
        ob.G = Gnufft_v2(args);
    else
        ob.G = Gnufft_v2(mask(:),args);
    end
elseif gridflag == 1
    ob.G = Gfft(N1,N2);
    if ~(K1 == N1)
        sprintf('Oversmapling not used with epiflag')
    end
% elseif gridflag == 2
%     if ~exist('mask','var')
%         ob.G = Gnufft_v2(args);
%     else
%         ob.G = Gnufft_v2(mask(:),args);
%     end
else %assume noncartesian in all dimensions
    if ~exist('mask','var')
        ob.G = Gnufft_v2(args);
    else
        ob.G = Gnufft_v2(mask(:),args);
    end
end

ob.is_masked = ob.G.is_masked;
ob.mask = ob.G.mask;

if (sum(abs(we(:))) == 0)
    sprintf('Setting L=1') 
    L=1;
end
% if (isvar('we_histo'))&(int_opt==2)
%     if (sum(abs(we_histo(:,1))) == 0)
%         sprintf('Setting L=0')
%         L=0;
%     end
%     ob.int = int_tim_seghanning(tt,L);
% elseif (isreal(we))
%     ob.int = int_tim_seghanning(tt,L);
% end
if (isvar('we_histo'))&&(int_opt==2)
    if (sum(abs(we_histo(:,1))) == 0)
        sprintf('Setting L=1')
        L=1;
    end
    ob.int = int_tim_seg(tt,L,we(:),int_opt,we_histo);
elseif (isreal(we))
    %ob.int = int_tim_seg(tt,L,we(:),int_opt);
    ob.int = int_tim_seghanning(tt,L);
elseif (~isreal(we))
    ob.int = int_tim_seghanning(tt,L);
end

ob.ttsegments = length(kx(:))/length(tt(:));

%This part added for T2*
if 0
    if ~(sum(abs(imag(we(:))))<eps)
        sprintf('Using R2* interpolator.')
        %keyboard
        ob.int = timr2_seg_int(tt,L,real(we(:)),-imag(we(:)),3, 1000);
    end
end

if flag_swt
    ob.flgswt = 1;
    ob.swt = sinc(kx/N).*sinc(ky/N);
else
    ob.flgswt = 0;
end

ob.we = we;
ob.tt = tt;
ob.is.empty	= logical(0);
%	ob.m = size(we);	% image size
%	ob.n = size(we);

ob = class(ob, 'fast_mr_v2');






