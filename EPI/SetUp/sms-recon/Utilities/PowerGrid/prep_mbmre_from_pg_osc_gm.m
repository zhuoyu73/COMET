function prep_mbmre_from_pg_osc

% hardwired to 16 slices; 24 reps
% need to make automatic

%for jj = 1:24
clear img
    for ii = 1:24
        if exist(sprintf('img_Slice%i_Avg0_Echo0_Phase0_mag.nii',ii-1))
        tmp1 = load_untouch_nii(sprintf('img_Slice%i_Avg0_Echo0_Phase0_mag.nii',ii-1));
        tmp2 = load_untouch_nii(sprintf('img_Slice%i_Avg0_Echo0_Phase0_phs.nii',ii-1));
        tmp3(:,:,:,:,ii) = tmp1.img.*exp(1i*tmp2.img);
        img = permute(tmp3,[1 2 3 5 4]);
        
else
        end
    end
          save img.mat img

    

mkdir recon
!mv img_* recon/