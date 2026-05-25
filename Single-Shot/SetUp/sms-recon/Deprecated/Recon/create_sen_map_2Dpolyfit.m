function sen = create_sen_map_2Dpolyfit(dirnm,recoInfo,flag_average)
% function sen = create_sen_map(dirnm)
%
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

if 1
   for coilIndex = 1:num_coils
      tmp = coil_imgs(:,:,:,coilIndex)./sos_img;
      tmp(find(isnan(tmp))) = 0;
      tmp(find(isinf(tmp))) = 0;
      coil_imgs(:,:,:,coilIndex) = tmp;
   end
end


mask = (abs(sos_img) > (0.1*max(abs(sos_img(:)))));

sos_img(find(isnan(sos_img))) = 0;

slice_mean = squeeze(mean(mean(abs(mask),2),1));
sl_mask = slice_mean>0;


for sliceIndex = 1:nz
    if sl_mask(sliceIndex)
       for coilIndex = 1:num_coils
               %sen(:,:,coilIndex) = sen(:,:,coilIndex)./abs(sos_img);
               sen(:,:,sliceIndex,coilIndex) = twodpolyfit(coil_imgs(:,:,sliceIndex,coilIndex), mask(:,:,sliceIndex), 3);
       end
    else
       for coilIndex = 1:num_coils
            sen(:,:,sliceIndex,coilIndex) = 1e4;
       end
    end
    
end

cd(curdir)

%reshape sen matrix for storage
save sen sen




