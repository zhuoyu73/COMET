function FM = create_field_map_telist(dirnm,recoInfo,TElist);


addpath('/shared/mrfil/code/MATLAB/Fessler_IRT/irt/')
setup_fessler_IRT

curdir = pwd;

% function FM = create_field_map(dirnm_FM1,dirnm_FM2,recoInfoFM1,recoInfoFM2);

nsl = recoInfo.nsl;

FM_tp_ind = 1:length(TElist);

% Load in images
for FM_tp_ind = 1:length(TElist)

  if ~exist(sprintf('%s/Recon/sen_%05d',dirnm,FM_tp_ind),'file')
     load(sprintf('%s/Recon/sen_%05d',dirnm,FM_tp_ind))
     img(find(isnan(img)))=0;
     eval(sprintf('Image%d = img;',FM_tp_ind))
  else
     sprintf('Did not find %s/Recon/sen_%05d \n',dirnm,FM_tp_ind)
     return
  end
end

slice_mean = squeeze(mean(mean(abs(Image1),2),1));
sl_mask = slice_mean>0;

 %FM = -mri_field_map_reg(cat(4,Image1(:,:,sl_mask), Image2(:,:,sl_mask)), [TE1 TE2],'l2b',-3);

data = cat(4,Image1(:,:,sl_mask), Image2(:,:,sl_mask));
for FM_tp_ind = 3:length(TElist)
    eval(sprintf('data = cat(4,data,Image%d(:,:,sl_mask));',FM_tp_ind))
end

 if (size(data,3) == 1)
     data = squeeze(data);
 end

 mask_img = abs(data(:,:,:,end))>(0.1*max(abs(col(data(:,:,:,end))))); %change by JLH 2013_09_05
 %mask_img = abs(Image1)>(0.1*max(abs(Image1(:))));
 FM = -mri_field_map_reg(data, TElist,'l2b',-3,'mask',mask_img);

 
 FMout = zeros(size(Image1));
 FMout(:,:,sl_mask) = FM;
 FM = FMout;
 
cd(curdir)
save FM FM

