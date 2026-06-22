function proc_mbmre_2

load t2stack.mat
load t2mask_bet.mat
load maskx.mat

mask = mask.*abs(1-maskx);
%mask = maskx;
save t2mask_final.mat mask

cd mre_process
save t2mask.mat mask

load mreimages.mat image
[ny,nx,nz,ni,nj] = size(image);
nk = ni*nj;

% note: will have to come up with a better way of making these automatic
dx = 240/nx;
dy = 240/ny;
dz = dx;
freq = 50;

image(isnan(image))=0;
imrs = reshape(image,[ny nx nz nk]);

mask_rep = repmat(mask,[1 1 1 nk]);
mask_nii = make_nii(mask_rep);
save_nii(mask_nii,'mre_mask.nii')

phs_nii = make_nii(mask_rep.*angle(imrs));
save_nii(phs_nii,'mre_phs.nii')

mag_rep = repmat(t2stack,[1 1 1 nk]);
mag_nii = make_nii(mag_rep);
save_nii(mag_nii,'mre_mag.nii')

!$FSLDIR/bin/prelude -a mre_mag.nii -p mre_phs.nii -o mre_output.nii -m mre_mask.nii
!gunzip -f mre_output.nii.gz

tmp = load_nii('mre_output.nii');
pimg = double(reshape(tmp.img,[ny nx nz ni nj]));
save mreimages_unwrap.mat pimg

%OSS_SNR = oss_snr_filter(pimg,[dx dy dz],freq,mask);

disp_img = flip(pimg,5)*1.764;
fft_img = fft(disp_img,[],5)/(nj/2);
wave_img = fft_img(:,:,:,:,2);

Zmotion = wave_img(:,:,:,1);
Ymotion = wave_img(:,:,:,2);
Xmotion = wave_img(:,:,:,3);

mreParams.subj = 'mbmre_for_inversion';
mreParams.FOVx = nx*dx;
mreParams.FOVy = ny*dy;
mreParams.FOVz = nz*dz;
mreParams.nx = ny;
mreParams.ny = ny;
mreParams.nz = nz;
mreParams.freq = freq;
%mreParams.oss_snr = OSS_SNR;

save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack')
cd ..
save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack')
