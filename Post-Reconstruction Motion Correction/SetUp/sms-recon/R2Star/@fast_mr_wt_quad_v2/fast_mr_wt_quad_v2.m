function ob = fast_mr_wt_quad_v2(kx,ky,kz,fov,Nx,Ny,Nz,Kx,Ky,Kz,J,tt,hanning_flag,we,r2,flag_swt,L,int_opt,we_histo,epiflag,nx,ny,nz,zop,gz,gridflag,mask)
%function ob = fast_mr_wt(kx,ky,kz,fov,Nx,Ny,Nz,Kx,Ky,Kz,J,tt,we,flag_swt,L,int_opt,we_histo,epiflag,nx,ny,nz,zop,gz)
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
ob.Nx = 0;
ob.Ny = 0;
ob.Nz = 0;
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

% if ((nargin < 18) || (nargin > 20))
%     help fast_mr_wt_quad_v2
%     error nargin
% end

%	fill object
Jm1 = J;
Jm2 = J;
if Nz == 1
    Jm3 = 1;
elseif Nz < 8
    Jm3 = 4;
else
    Jm3 = J;
end

if Nz == 1
    if nz ~= 1
        error('Cannot use oversampling in z if a 2D transform is used.');
    end
end


% % ob.st = nufft2_init([2*pi/(Nx)*ky,2*pi/(Ny)*kx],Nx,Ny,Jm1,Jm2,Kx,Ky,[N/2-1/2,N/2-1/2]*(1),[],'kb');
% %ob.st = nufft_init([2*pi/(Nx)*ky,2*pi/(Ny)*kx],[Nx,Ny],[Jm1,Jm2],[Kx,Ky],[Nx/2-1/2,Ny/2-1/2]*(1));
% args = {[2*pi/Nx*ky,2*pi/Ny*kx],[Nx,Ny],[Jm1,Jm2],[Kx,Ky],[Nx/2-1/2,Nx/2-1/2]*(1)};
% %args = {[2*pi/Nx*ky,2*pi/Ny*kx],[Nx,Ny],[Jm1,Jm2],[Kx,Ky],[N/2,N/2]*(1)};
%
% if ~epiflag
%     ob.G = Gnufft(args);
% else
%     ob.G = Gfft(Nx/nx,Nx);
%     if ~(K == N)
%         sprintf('Oversmapling not used with epiflag')
%     end
% end

if(length(kz) > 1)
    warning('Performing 3D NUFFT...')
    if gridflag == 2
        args = {[2*pi/(Ny*ny)*ky,2*pi/(Nx*nx)*kx, 2*pi/(Nz*nz)*kz],[Nx*nx,Ny*ny, Nz*nz],[Jm1,Jm2,Jm3],[nx*Kx,ny*Ky,nz*Kz],[nx*Nx/2-1/2,ny*Ny/2-1/2, nz*Nz/2-1/2], 'table', 2^10, 'minmax:kb'};
        warning('Using table lookup...');
    else
        args = {[2*pi/(Ny*ny)*ky,2*pi/(Nx*nx)*kx, 2*pi/(Nz*nz)*kz],[nx*Nx,nx*Ny, nz*Nz],[Jm1,Jm2,Jm3],[nx*Kx,ny*Ky,nz*Kz],[nx*Nx/2-1/2,ny*Ny/2-1/2, nz*Nz/2-1/2], 'minmax:kb', 2^10};
    end
    %     args = [ {[2*pi/Nx*ky,2*pi/Ny*kx, 2*pi/Nz*kz],[Nx,Ny, Nz],[Jm1,Jm2,Jm3],[Kx,Ky,Kz]}, args];
    %     args = {[2*pi/Nx*ky,2*pi/Ny*kx, 2*pi/Nz*kz],[Nx,Ny, Nz],[Jm1,Jm2,Jm3],[Kx,Ky,Kz],[Nx/2,Ny/2, Nz/2], 'table', 2^10, 'minmax:kb'};
    %     args = {[2*pi/Nx*ky,2*pi/Ny*kx, 2*pi/Nz*kz],[Nx,Ny, Nz],[Jm1,Jm2,Jm3],[Kx,Ky,Kz],[Nx/2,Ny/2, Nz/2], 'minmax:kb', 2^10};
    %     args = {[2*pi/Nx*ky,2*pi/Ny*kx, 2*pi/Nz*kz],[Nx,Ny, Nz],[Jm1,Jm2,Jm3],[Kx,Ky,Kz],[Nx/2-1/2,Ny/2-1/2, Nz/2-1/2], 'linear', 2^10};
else
    if gridflag == 2
        args = {[2*pi/Nx*ky,2*pi/Ny*kx],[nx*Nx,ny*Ny],[Jm1,Jm2],[nx*Kx,ny*Ky],[Nx/2-1/2,Ny/2-1/2]*(1), 'table', 2^10, 'minmax:kb'};
        warning('Using table lookup...');
    else
        args = {[2*pi/Nx*ky,2*pi/Ny*kx],[nx*Nx,ny*Ny],[Jm1,Jm2],[nx*Kx,ny*Ky],[Nx/2-1/2,Ny/2-1/2]*(1), 'minmax:kb', 2^10};
    end
    
end

if (gridflag == 0)||(gridflag == 2)
    if ~exist('mask','var')
        ob.G = Gnufft_v2(args);
    else
        ob.G = Gnufft_v2(mask(:),args);
    end
elseif gridflag == 1
    ob.G = Gfft(Nx,Ny);
    if ~(Kx == Nx)
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


if (ndims(we) == 2)
    we = reshape(we,Nx,Ny,Nz); % issue with jointest and four_echo_mri, field map is in column
    r2 = reshape(r2,Nx,Ny,Nz);
end


H1 = sparse(kron(eye(Nx),1/(ny)*ones(ny,1)));
H2 = sparse(kron(eye(Ny),1/(nx)*ones(nx,1)));
H3 = sparse(kron(eye(Nz),1/(nz)*ones(nz,1)));
HH = sparse(kron(H3,kron(H2,H1)));

%Resize Field map
we_out = zeros(ny*Ny,nx*Nx,Nz*nz);
r2_out = zeros(ny*Ny,nx*Nx,Nz*nz);


% if (ndims(we) == 2)
%     nzwe = 1;
%     szx = size(we,2);
%     szy = size(we,1);
% else
nzwe = size(we,3);
szx = size(we,2);
szy = size(we,1);
% end
wem = zeros(szy,szx,nzwe);

%Interpolate field map
if Nz == 1
    if ~(ndims(we) == 2)
        error('Field map does not match transform size in z');
    else
        wec = we;
    end
    
    [nwey,nwex] = size(wec);
    
    if ~((nwey == ny*Ny)&&(nwex == nx*Nx))
        [wxl1,wyl1]= meshgrid([-(nwex)/2:(nwex)/2-1]./(nwex),[-(nwey)/2:(nwey)/2-1]./(nwey));
        [wxl2,wyl2]= meshgrid([-(nx*Nx)/2:(nx*Nx)/2-1]./(nx*xN),[-(ny*Ny)/2:(ny*Ny)/2-1]./(ny*Ny));
        wen = interp2(wxl1+1/(2*nwex),wyl1+1/(2*nwey),wec,wxl2+(1/(2*nx*N)),wyl2+(1/(2*ny*N)),'cubic');
        r2n = interp2(wxl1+1/(2*nwex),wyl1+1/(2*nwey),r2,wxl2+(1/(2*nx*N)),wyl2+(1/(2*ny*N)),'cubic');
        %       r2n = reshape(HH*r2(:),nx*N,ny*N);
    else
        wen = wec;
        r2n=r2;
    end
    
    we_out(:,:) = wen;
    r2_out(:,:) = r2n;
    
else % 3D Case
    [nwey, nwex, nwez] = size(we);
    
    if ~(((nwey == ny*Ny)&&(nwex == nx*Nx))&&(nwez == nz*Nz))
        [wxl1,wyl1,wzl1]= meshgrid([-(nwex)/2:(nwex)/2-1]./(nwex),[-(nwey)/2:(nwey)/2-1]./(nwey),[-(nwez)/2:(nwez)/2-1]./(nwez));
        [wxl2,wyl2,wzl2]= meshgrid([-(nx*Nx)/2:(nx*Nx)/2-1]./(nx*Nx),[-(ny*Ny)/2:(ny*Ny)/2-1]./(ny*Ny),[-(nz*Nz)/2:(nz*Nz)/2-1]./(nz*Nz));
        wen = interp3(wxl1+1/(2*nwex),wyl1+1/(2*nwey),wzl1+1/(2*nwez),we,wxl2+(1/(2*nx*Nx)),wyl2+(1/(2*ny*Ny)),wzl2+(1/(2*nz*Nz)),'cubic');
        r2n = interp3(wxl1+1/(2*nwex),wyl1+1/(2*nwey),wzl1+1/(2*nwez),r2,wxl2+(1/(2*nx*Nx)),wyl2+(1/(2*ny*Ny)),wzl2+(1/(2*nz*Nz)),'cubic');
    else
        wen = we;
        r2n = r2;
    end
    we_out = wen;
    r2_out = r2n;
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
%ob.N = N;

if flag_swt
    ob.flgswt = 1;
    if Nz > 1
        ob.swt = sinc(kx/Nx).*sinc(ky/Ny).*sinc(kz/Nz);
    else
        ob.swt = sinc(kx/Nx).*sinc(ky/Ny);
    end
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

ob = class(ob, 'fast_mr_wt_quad_v2');
