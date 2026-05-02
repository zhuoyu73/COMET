function ob = fast_mr(kx,ky,fov,N,K,J,tt,we,flag_swt,L,int_opt,we_histo, epiflag)
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
%    epiflag = this is 1 when EPI sequence.  No need to oversample, just fft.  
% Note:  if we = zeros, then L=0 is used.
%  Brad Sutton   Jeff Fessler 
%  University of Michigan
%  June 17, 2002 
%
% History:  added epiflag 4/9/04
  
% This is an attempt to make MRI object valid for sensitivity
%   encoded runs also
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
%ob.version = 1.0;

if nargin == 0
	ob = class(ob, 'fast_mr');
	return
end

if isa(kx, 'fast_mr')
	ob = kx;
        if nargin == 1
           return
        elseif (isstr(ky))
           eval(sprintf('ob.%s = fov;',ky)) 
        end
	return
end

if (nargin < 11)|(nargin>13)
	help fast_mr
	error nargin
end

if ~exist('epiflag','var')
  epiflag = 0;
end

	%	fill object
        N1 = N;
        N2 = N;
        Jm1 = J;
        Jm2 = J;
        K1 = K; 
        K2 = K;
        %ob.st =
        %nufft2_init([2*pi/N1*ky,2*pi/N2*kx],N1,N2,Jm1,Jm2,K1,K2,[N/2-1/2,N/2-1/2]*(1),[],'kb');
	args = {[2*pi/N1*ky,2*pi/N2*kx],[N1,N2],[Jm1,Jm2],[K1,K2],[N/2-1/2,N/2-1/2]*(1)};
	%args = {[2*pi/N1*ky,2*pi/N2*kx],[N1,N2],[Jm1,Jm2],[K1,K2],[N/2,N/2]*(1)};
	if ~epiflag
          ob.G = Gnufft(args);
        else
          ob.G = Gfft(N,N);
          if ~(K == N)
            sprintf('Oversmapling not used with epiflag')
          end
        end
        
	if (sum(abs(we(:))) == 0)
             sprintf('Setting L=0') 
             L=0;
         end

         if (isvar('we_histo'))&(int_opt==2)
            if (sum(abs(we_histo(:,1))) == 0)
                sprintf('Setting L=0')
                L=0;
            end

            ob.int = int_tim_seg(tt,L,we(:),int_opt,we_histo);
         elseif (isreal(we))
            ob.int = int_tim_seg(tt,L,we(:),int_opt);
         end

 %This part added for T2*
        if ~(sum(abs(imag(we(:))))<eps)
	    sprintf('Using R2* interpolator.')
	    %keyboard
            ob.int = timr2_seg_int(tt,L,real(we(:)),-imag(we(:)),3, 1000);
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

	ob = class(ob, 'fast_mr');






