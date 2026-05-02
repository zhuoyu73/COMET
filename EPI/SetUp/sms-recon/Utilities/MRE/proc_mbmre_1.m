function proc_mbmre_1

load img.mat

imgx = squeeze(img);
for ii = 1:size(img,4)
    imgy(:,:,ii+(0:size(img,4):3.*size(img,4)),:) = squeeze(imgx(:,:,:,ii,:));
end

imgy = flip(imgy,1);

[ny, nx, nz, nk] = size(imgy);
dy = 240/ny;
dx = 240/nx;
dz = dy;
freq = 50;

imgraw = reshape(imgy,[ny nx nz 2 3 nk/6]);

save imgraw_mbmre.mat imgraw

t2stack = mean(mean(mean(abs(imgraw),6),5),4);
t2nii = make_nii(flipdim(flipdim(permute(t2stack,[2 1 3]),1),2),[dy dx dz]);
save_nii(t2nii,'t2stack.nii')

!$FSLDIR/bin/bet2 t2stack.nii t2bet.nii -m -v -f 0.4 -w 1
!gunzip -f t2bet_mask.nii.gz
!gunzip -f t2bet.nii.gz
!cp t2bet_mask.nii t2mask.nii
!rm t2bet_mask.nii

tmp = load_nii('t2mask.nii');
seD1 = strel('diamond',2);
mask = imerode(double(permute(flipdim(flipdim(tmp.img,1),2),[2 1 3])),seD1);
save t2mask_bet.mat mask
save t2stack.mat t2stack

mkdir mre_process
cd mre_process

image(:,:,:,1,:) = squeeze(imgraw(:,:,:,1,3,:)./imgraw(:,:,:,2,3,:)); % Z
image(:,:,:,2,:) = squeeze(imgraw(:,:,:,1,1,:)./imgraw(:,:,:,2,1,:)); % Y
image(:,:,:,3,:) = squeeze(imgraw(:,:,:,1,2,:)./imgraw(:,:,:,2,2,:)); % X

save mreimages.mat image
cd ..

