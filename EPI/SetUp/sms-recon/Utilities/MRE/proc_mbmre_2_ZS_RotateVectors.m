
load('affine_matrices_crosscor_ZS.mat');
motions_tmp = squeeze(affine_matrices(:,:,1,:,:,1,1,1,1,:));
motions = squeeze(mean(mean(motions_tmp, 3), 4));
affines = reshape(motions, [3 3 2 3 4]);
affines2 = permute(affines, [1 2 5 3 4]);

load t2stack.mat
load t2mask_bet.mat
load maskx.mat
mask = mask.*abs(1-maskx);
save t2mask_final.mat mask

load imgraw_mbmre.mat imgraw 

mkdir mre_process_img_reg_wave
cd mre_process_img_reg_wave
mkdir dcfiles

clear imgraw2
clear pimg
imgraw_tmp = permute(imgraw,[1 2 3 6 4 5]);

imgraw2(:,:,:,:,:,1)=imgraw_tmp(:,:,:,:,:,3);
imgraw2(:,:,:,:,:,2)=imgraw_tmp(:,:,:,:,:,1);
imgraw2(:,:,:,:,:,3)=imgraw_tmp(:,:,:,:,:,2);

[ny,nx,nz,ni,nj,nl] = size(imgraw2);
nk = ni*nj;

dx = 240/nx;
dy = 240/ny;
dz = dx;
freq = 50;


phaseOffsets = [0 pi/2 pi 3*pi/2];
polSigns = [+1 -1];
E_nom = eye(3);

[Nx,Ny,Nz,Nphase,Npol,Ndir] = size(imgraw2);
U_phase = zeros(Nx,Ny,Nz,3,Nphase);

for iph = 1:Nphase
    m_all = zeros(Nx,Ny,Nz,Ndir);
    E_all = zeros(Ndir,3);

    for idir = 1:Ndir
        e_nom = E_nom(idir,:).'; 
        phase_pol = zeros(Nx,Ny,Nz,Npol);
        E_pol = zeros(Npol,3);

        for ipol = 1:Npol
            phase_pol(:,:,:,ipol) = angle(imgraw2(:,:,:,iph,ipol,idir));
            M = affines2(:,:,iph,ipol,idir);
            A = M(1:2,1:2);

            [U_svd,~,V_svd] = svd(A);
            R2 = U_svd * V_svd';

            if det(R2) < 0
                U_svd(:,end) = -U_svd(:,end);
                R2 = U_svd * V_svd';
            end
            R = eye(3);
            R(1:2,1:2) = R2;
            e_eff = R' * e_nom;
            E_pol(ipol,:) = e_eff.';
        end


        phi_diff = squeeze(imgraw2(:,:,:,iph,1,idir)./imgraw2(:,:,:,iph,2,idir));
        E_eff_dir = mean(E_pol,1);

        m_all(:,:,:,idir) = phi_diff;
        E_all(idir,:) = E_eff_dir;
    end

    Einv = pinv(E_all);

    mx = m_all(:,:,:,1);
    my = m_all(:,:,:,2);
    mz = m_all(:,:,:,3);

    U_phase(:,:,:,1,iph) = Einv(1,1)*mx + Einv(1,2)*my + Einv(1,3)*mz;
    U_phase(:,:,:,2,iph) = Einv(2,1)*mx + Einv(2,2)*my + Einv(2,3)*mz;
    U_phase(:,:,:,3,iph) = Einv(3,1)*mx + Einv(3,2)*my + Einv(3,3)*mz;

    save mreimagesUphase.mat U_phase

    [ny,nx,nz,ni] = size(squeeze(U_phase(:,:,:,:,iph)));
    nk=ni;
    U_phase2 = U_phase(:,:,:,:,iph);
    U_phase2(isnan(U_phase2))=0;
    imrs = reshape(U_phase2,[ny nx nz nk]);
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

mreParams.subj = 'mbmre_for_inversion';
mreParams.FOVx = nx*dx;
mreParams.FOVy = ny*dy;
mreParams.FOVz = nz*dz;
mreParams.nx = ny;
mreParams.ny = ny;
mreParams.nz = nz;
mreParams.freq = freq;
mreParams.oss_snr = OSS_SNR;

cd ../

save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack','OSS_SNR')
