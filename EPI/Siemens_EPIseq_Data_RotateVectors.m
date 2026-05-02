%% Siemens_EPIseq_Data_RotateVectors.m 
% Processing Part 1
mag_nii = load_untouch_nii('mag.nii');
phs_nii = load_untouch_nii('phs.nii');
magimg = mag_nii.img;
phsimg = phs_nii.img;

mag = permute((flip(magimg,2)),[2 1 3 4]);
phs = permute((flip(phsimg,2)),[2 1 3 4]);

dirlist_mag = dir('mag/*.nii*');
dirlist_phs = dir('phs/*.nii*');

[ny,nx,nz,~] = size(mag);
dx = mag_nii.hdr.dime.pixdim(2); 
dy = mag_nii.hdr.dime.pixdim(3);
dz = mag_nii.hdr.dime.pixdim(4);

dirs = 3; %number of directions (y, x, z) 
posneg = 2; %positive or negative polarity 
offset = 4; %phase offsets 
freq = 50;  

mag = double(mag);
phs = pi*((double(phs)-2048)/2048); 
cplx_img = mag.*exp(1i*phs);

imgraw = reshape(cplx_img,[ny nx nz posneg dirs offset]); 

save imgraw_ep2d.mat imgraw

t2stack = mean(mean(mean(abs(imgraw),6),5),4);  


t2nii = make_nii(flipdim(flipdim(permute(t2stack,[2 1 3]),1),2),[dy dx dz]);  
save_nii(t2nii,'t2stack.nii')


!$FSLDIR/bin/bet2 t2stack.nii t2bet.nii -m -v -f 0.25 -w 1.3
!gunzip -f t2bet.nii_mask.nii.gz
!gunzip -f t2bet.nii.gz
!cp t2bet.nii_mask.nii t2mask.nii
!rm t2bet.nii_mask.nii

tmp = load_nii('t2mask.nii');
seD1 = strel('diamond',2); 
mask = imerode(double(permute(flipdim(flipdim(tmp.img,1),2),[2 1 3 4])),seD1);
save t2mask_bet.mat mask
save t2stack.mat t2stack

%% Manual Masking

load('t2mask_bet.mat')
load('t2stack.mat')
tmp = (t2stack.*mask)./(1000);

maskx = zeros(size(mask));

% look at the whole brain
figure;im(tmp);caxis([0 1])

% look at the mask you have created (the negative mask)
figure;im(tmp.*mask.*abs(1-maskx));caxis([0 1])


for ss = 48:-1:1
    ss
    maskx(:,:,ss) = double(roipoly(tmp(:,:,ss)));
end
save maskx.mat maskx 

%% Process Part 2 
clear image

load maskx.mat

mask = mask.*abs(1-maskx);

save t2mask_final.mat mask

%% Motion Correction
load('imgraw_ep2d.mat') 
SubjectName = 'imgraw_ep2d'; 
tp = 24;
mkdir('individual_nii')
imgraw2 = reshape(imgraw,[size(imgraw,1) size(imgraw,2) size(imgraw,3) tp]); 
for ii=1:size(imgraw2,4)
    imgraw_rep1 = imgraw2(:,:,:,ii);
    tmp = load_nii('t2stack.nii');
    tmp.img= imgraw_rep1;
    cd('individual_nii')
    save_nii(tmp,sprintf('rep%d.nii',ii))
    cd ..
end
%imgraw mag
for ii=1:size(imgraw2,4)
    imgraw_rep1 = imgraw2(:,:,:,ii);
    tmp = load_nii('t2stack.nii');
    tmp.img= abs(imgraw_rep1);
    cd('individual_nii')
    save_nii(tmp,sprintf('rep%d_mag.nii',ii))
    cd ..
end
%imgraw real 
for ii=1:size(imgraw2,4)
    imgraw_rep1 = imgraw2(:,:,:,ii);
    tmp = load_nii('t2stack.nii');
    tmp.img= real(imgraw_rep1);
    cd('individual_nii')
    save_nii(tmp,sprintf('rep%d_real.nii',ii))
    cd ..
end
%imgraw imag
for ii=1:size(imgraw2,4)
    imgraw_rep1 = imgraw2(:,:,:,ii);
    tmp = load_nii('t2stack.nii');
    tmp.img= imag(imgraw_rep1);
    cd('individual_nii')
    save_nii(tmp,sprintf('rep%d_imag.nii',ii))
    cd ..
end

    cd('individual_nii')

    
for ii=1:tp
eval(sprintf('!$FSLDIR/bin/flirt -in rep%d_mag.nii -ref rep1_mag.nii -out rep%d_reg_mag -omat rep%d_reg_mag.mat -dof 6',ii,ii,ii))
end
!gunzip -f *.nii.gz

for ii=1:tp
eval(sprintf('!$FSLDIR/bin/flirt -in rep%d_real.nii -ref rep1_reg_mag.nii -out rep%d_reg_real.nii -init rep%d_reg_mag.mat -applyxfm',ii,ii,ii))
eval(sprintf('!$FSLDIR/bin/flirt -in rep%d_imag.nii -ref rep1_reg_mag.nii -out rep%d_reg_imag.nii -init rep%d_reg_mag.mat -applyxfm',ii,ii,ii))
end
!gunzip -f *.nii.gz


for ii=1:tp
    tmp_real = load_nii(sprintf('rep%d_reg_real.nii',ii));
    tmp_imag = load_nii(sprintf('rep%d_reg_imag.nii',ii));
    imgraw_new(:,:,:,ii) = (tmp_real.img+(1i*tmp_imag.img));
end

imgraw = reshape(imgraw_new,[size(imgraw,1) size(imgraw,2) size(imgraw,3) 2 3 4]); 
cd ..
save imgraw_ep2d_MoCo.mat imgraw

%% Wave Filed Correction
load('imgraw_ep2d_MoCo.mat')
load t2mask_final.mat
load t2stack.mat
cd ('individual_nii/')
motions = zeros([4 4 24]);
for ii=1:24
motions(:,:,ii) = load(sprintf('rep%d_reg_mag.mat',ii),'-ascii');
end

affines = reshape(motions,[4 4 2 3 4]);
affines2 = permute(affines,[1 2 5 3 4]);

cd ..

phaseOffsets = [0 pi/2 pi 3*pi/2];
E_nom = eye(3); 

imgraw_tmp = permute(imgraw,[1 2 3 6 4 5]);
clear imgraw2
imgraw2(:,:,:,:,:,1)=imgraw_tmp(:,:,:,:,:,3);
imgraw2(:,:,:,:,:,2)=imgraw_tmp(:,:,:,:,:,1);
imgraw2(:,:,:,:,:,3)=imgraw_tmp(:,:,:,:,:,2);

[Nx,Ny,Nz,Nphase,Npol,Ndir] = size(imgraw2);
U_phase = zeros(Nx,Ny,Nz,3,Nphase);


for iph = 1:Nphase
    m_all = zeros(Nx,Ny,Nz,Ndir);
    E_all = zeros(Ndir,3);

    for idir = 1:Ndir
        e_nom = E_nom(idir,:).';
        E_pol = zeros(Npol,3);

        for ipol = 1:Npol
            M = affines2(:,:,iph,ipol,idir); 
            A = M(1:3,1:3);

                    [U_svd,~,V_svd] = svd(A);
                    R = U_svd * V_svd';
        
                    if det(R) < 0
                        U_svd(:,end) = -U_svd(:,end);
                        R = U_svd * V_svd';
                    end

            e_eff = R' * e_nom;
            E_pol(ipol,:) = e_eff.';
        end

        phi_diff = squeeze((imgraw2(:,:,:,iph,1,idir))./(imgraw2(:,:,:,iph,2,idir))); 
        E_eff_dir = mean(E_pol,1); 

        m_all(:,:,:,idir) = phi_diff;
        E_all(idir,:) = E_eff_dir;
    end
    
    disp(E_all)
    disp(det(E_all))
    disp(cond(E_all))

    Einv = pinv(E_all);  

    mx = m_all(:,:,:,1);
    my = m_all(:,:,:,2);
    mz = m_all(:,:,:,3);

    U_phase(:,:,:,1,iph) = Einv(1,1)*mx + Einv(1,2)*my + Einv(1,3)*mz;
    U_phase(:,:,:,2,iph) = Einv(2,1)*mx + Einv(2,2)*my + Einv(2,3)*mz;
    U_phase(:,:,:,3,iph) = Einv(3,1)*mx + Einv(3,2)*my + Einv(3,3)*mz;
    save mreimagesUphase.mat U_phase

[nx,ny,nz,ni] = size(squeeze(U_phase(:,:,:,:,iph)));
nk=ni;
nj=4;
U_phase2 = U_phase(:,:,:,:,iph);
U_phase2(isnan(U_phase2))=0;
imrs = reshape(U_phase2,[nx ny nz nk]);
save t2mask.mat mask

mask_rep = repmat(mask,[1 1 1 nk]);
mask_nii = make_nii(mask_rep);
save_nii(mask_nii,'mre_mask.nii')

phs_nii = make_nii(mask_rep.*angle(imrs));
save_nii(phs_nii,sprintf('mre_phs%d.nii',iph))

mag_rep = repmat(t2stack,[1 1 1 nk]);
mag_nii = make_nii(mag_rep);
save_nii(mag_nii,sprintf('mre_mag%d.nii',iph))

tic
eval(sprintf('!$FSLDIR/bin/prelude -a mre_mag%d.nii -p mre_phs%d.nii -o mre_output%d.nii -m mre_mask.nii -v',iph,iph,iph))
!gunzip -f *.nii.gz
toc

tmp = load_nii(sprintf('mre_output%d.nii',iph));
pimg(:,:,:,:,iph) = double(tmp.img);
end

save mreimages_unwrap.mat pimg

OSS_SNR = oss_snr_filter(pimg,[dx dy dz],freq,mask);

disp_img = flip(pimg,5)*1.764;
fft_img = fft(disp_img,[],5)/(nj/2);
wave_img = fft_img(:,:,:,:,2);

Zmotion = wave_img(:,:,:,1);
Ymotion = wave_img(:,:,:,2);
Xmotion = wave_img(:,:,:,3);

mreParams.subj = 'mre_for_inversion';
mreParams.FOVx = nx*dx;
mreParams.FOVy = ny*dy;
mreParams.FOVz = nz*dz;
mreParams.nx = nx;
mreParams.ny = ny;
mreParams.nz = nz;
mreParams.freq = freq;
mreParams.oss_snr = OSS_SNR;

save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack','OSS_SNR')
