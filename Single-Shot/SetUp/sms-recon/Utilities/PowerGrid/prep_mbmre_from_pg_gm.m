function mergePowerGridFileOutput(NSlices,NReps)

%Converts PowerGrid Output to one Image

for ii = 1:NReps
    for jj = 1:NSlices
            tmp1 = load_untouch_nii(sprintf('img_Slice%i_Rep%i_Avg0_Echo0_Phase0_mag.nii',jj-1,ii-1),1,1,1,1);
            tmp2 = load_untouch_nii(sprintf('img_Slice%i_Rep%i_Avg0_Echo0_Phase0_phs.nii',jj-1,ii-1),1,1,1,1);
            img(:,:,:,jj,ii) = tmp1.img.*exp(1i*tmp2.img);
    end
end

save img.mat img

mkdir recon
!mv img_* recon/
