function FM = create_field_map_l2b0(dirnm_FM1,dirnm_FM2,recoInfoFM1,recoInfoFM2);

addpath('~/MRFIL/code/MATLAB/Fessler_IRT/irt/')
setup_fessler_IRT

curdir = pwd;

% function FM = create_field_map(dirnm_FM1,dirnm_FM2,recoInfoFM1,recoInfoFM2);

nsl = recoInfoFM1.nsl;
TE1 = recoInfoFM1.TE;
TE2 = recoInfoFM2.TE;

FM_tp_ind = 1;

% Load in images
if ~exist(sprintf('%s/Recon/sen_%05d',dirnm_FM1,FM_tp_ind),'file')
   load(sprintf('%s/Recon/sen_%05d',dirnm_FM1,FM_tp_ind))
   Image1 = img;
else
   sprintf('Did not find %s/Recon/sen_%05d \n',dirnm_FM1,FM_tp_ind)
   return
end

if ~exist(sprintf('%s/Recon/sen_%05d',dirnm_FM2,FM_tp_ind),'file')
   load(sprintf('%s/Recon/sen_%05d',dirnm_FM2,FM_tp_ind))
   Image2 = img;
else
   sprintf('Did not find %s/Recon/sen_%05d \n',dirnm_FM2,FM_tp_ind)
   return
end

Image1(find(isnan(Image1))) = 0;
Image2(find(isnan(Image2))) = 0;

slice_mean = squeeze(mean(mean(abs(Image1),2),1));
sl_mask = slice_mean>0;

 FM = -mri_field_map_reg(cat(4,Image1(:,:,sl_mask), Image2(:,:,sl_mask)), [TE1 TE2],'l2b',0);

 
 FMout = zeros(size(Image1));
 FMout(:,:,sl_mask) = FM;
 FM = FMout;
 
cd(curdir)
save FM FM

