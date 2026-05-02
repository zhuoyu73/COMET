function ob = LRobj_v2(A, v)
%function ob = fast_mr(kx,ky,fov,N,K,J,tt,we,flag_swt,L,int_opt,we_histo,epiflag)
%	Construct MRI object, which can do Ax and A'y operations
% typical input params:
%   kx goes form -Nx/2 to Nx/2:  unitless version   [nro nshots]
%   ky goes from -Ny/2 to Ny/2:  unitless version   [nro nshots]
%   kz goes from -Nz/2 to Nz/2:  unitless version   [nro nshots]
%    N = 64,128, MTX size
%    K oversampling   = 2*N, 3*N, 1.5*N  
%    J = number of neighbors, typ. 5,6,7,8,...
%    L number of time segments.   L=5,6,7,...
%    ntp is the number of images that will be used
%    we_histo is a histogram of field map (we) for int_opt=2
%            has two columns:
%                    column 1:  bin centers for histogram
%                    column 2:  histogram values at those bin centers
%    int_opt = 1 for using input field map to find min-max int.
%    int_opt = 2 for histogram and include we_histo for histogram of
% 	mask = logical mask of pixels 
%    field map 
% Note:  if we = zeros, then L=0 is used.
%
%  tt: timing vector for ONE shot (in a multi-shot) acquisition
% 
% 
%  Brad Sutton   Jeff Fessler 
%  University of Michigan
%  June 17, 2002 
%
% Joe Holtrop
% University of Illinois
% October 1, 2013
  


%	default object
ob.A = 0;
ob.is.empty	= logical(1);
ob.is.transpose = logical(0);

if nargin == 0
	ob = class(ob, 'LRobj_v2');
	return
end

%check vector sizes
if size(v,2) ~= length(A)
    keyboard
end

%	fill object
ob.A = A;
ob.v = v;
ob.ntp = size(v,2);
ob.R = size(v,1); %Rank of model being used
ob.is.empty	= logical(0);

ob = class(ob, 'LRobj_v2');






