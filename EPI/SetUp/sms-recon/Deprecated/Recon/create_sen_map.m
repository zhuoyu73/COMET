function sen = create_sen_map(dirnm,recoInfo,flag_average)
% function sen = create_sen_map(dirnm)
%
% Inputs
%    dirnm is the directory name where the coil images have been reconstructed
%

curdir = pwd;

if ~exist('flag_average','var')
   flag_average = 0;
end

if flag_average
   sprintf('Temporal averaging for sensitivity maps not implemented yet \n')
   sprintf('Type return to continue \n')
   keyboard
end


% First verify that the coil images are in the director, check dirnm and dirnm/Recon
cd(dirnm)
coil_list = dir('c01*.mat');
if isempty(coil_list)
   cd('Recon')
   coil_list = dir('c01*.mat');
end

if isempty(coil_list)
   sprintf('Coil images not found in %s nor in %s/Recon \n',dirnm,dirnm)
    keyboard
   cd(curdir)
   return
end


nx = recoInfo.N;
ny = recoInfo.N;
nz = recoInfo.nsl;
num_coils = recoInfo.num_coils;

coil_imgs = zeros(ny,nx,nz,num_coils);
sen = zeros(ny,nx,nz,num_coils);

sen_tp_index = 1;

for coilIndex = 1:num_coils
    load(sprintf('c%02d_%05d',coilIndex,sen_tp_index))
    coil_imgs(:,:,:,coilIndex) = img;
end

sos_img = sos_image(coil_imgs);

%create a mask
mask = (abs(sos_img) > (0.1*max(abs(sos_img(:)))));
[x y] = meshgrid(1:nx,1:ny);
mask_circ = sqrt((x-(nx+1)/2).^2+(y-(ny+1)/2).^2)<((nx)/2);
mask_circ = repmat(mask_circ,[1 1 nz]);
mask = mask.*mask_circ;

%mask coil images with a circular FOV
for coilIndex= 1:num_coils
    coil_imgs(:,:,:,coilIndex) = coil_imgs(:,:,:,coilIndex).*mask_circ;
end

%% SENSE map based on iterative smoothing
beta = 1e-6;
niter = 40;
init_img = 1e-2*ones(nx*ny,1);
Cn = Roughness_Penalty(ones(nx,ny))*beta;
for sliceIndex = 1:nz
    A = sensemult(mask(:,:,sliceIndex).*sos_img(:,:,sliceIndex));
    for coilIndex = 1:num_coils
        tmp = qpwls_pcg(init_img, A, 1,col(coil_imgs(:,:,sliceIndex,coilIndex)), 0, Cn, 1, niter);
        sen(:,:,sliceIndex,coilIndex) = reshape(tmp(:,niter),ny,nx);
    end
end


cd(curdir)

%reshape sen matrix for storage
save sen sen
end

function X = Roughness_Penalty(mask)
%function X = C3D_sparse(mask)
% Creates a roughness penalty matrix for the mask image.
% mask contains the penalty for each point in the image
% X is a sparse matrix containing penalties
        
    [szx szy szz] = size(mask);
    
    C_tmp = ones(szx,szy,szz);
    C_tmp2 =spdiag(col(C_tmp));
    mask1= spdiag(col(mask));
    
    C1 =C_tmp;
    C2 =C_tmp;
    C3 =C_tmp;
    
    %Set edges with no neighbors equal to zero
    C1(szx,:,:) = 0;
    C2(:,szy,:) = 0;
    C3(:,:,szz) = 0;
    
    %shift matrices to specify a neighbor
    C4 =circshift(C_tmp2,[-1 0]);
    C5 =circshift(C_tmp2,[-szy 0]);
    C6 =circshift(C_tmp2,[-szy*szx 0]);
    
    %C-tmp2 - C4 gives a C matrix with neighbors for all points
    %multiplying by C1 removes neighbors on the edges that do not have
    %neighbors
    %mask1  determines the weighting for each point
    C7=mask1*spdiag(col(C1))*(C_tmp2-C4);
    C8=mask1*spdiag(col(C2))*(C_tmp2-C5);
    C9=mask1*spdiag(col(C3))*(C_tmp2-C6);

    %concatanate all dimensions together
    C10=cat(1,C7,C8);
    X=cat(1,C10,C9);
end

