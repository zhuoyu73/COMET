function [FM,FMsmoothed] = estimFieldMap_alt(rInfo,FMImages,mask)
%ESTIMFIELDMAP_ALT Summary of this function goes here
%   Detailed explanation goes here

% se = strel('diamond',2);
% maskX = imerode(mask,se);
maskX = mask;

FMImages_mean = mean(abs(FMImages),4);
tmp = make_nii(FMImages_mean,[rInfo.FOV*10/rInfo.N rInfo.FOV*10/rInfo.N rInfo.sliceThickness]);
save_nii(tmp,'FMmag_tmp.nii')

 FMImages_sub = angle(FMImages(:,:,:,2)./FMImages(:,:,:,1));
%FMImages_sub = angle(FMImages(:,:,:,3)./FMImages(:,:,:,2));
tmp = make_nii(FMImages_sub,[rInfo.FOV*10/rInfo.N rInfo.FOV*10/rInfo.N rInfo.sliceThickness]);
save_nii(tmp,'FMphs_tmp.nii')

tmp = make_nii(double(maskX),[rInfo.FOV*10/rInfo.N rInfo.FOV*10/rInfo.N rInfo.sliceThickness]);
save_nii(tmp,'FMmsk_tmp.nii')

!$FSLDIR/bin/prelude -p FMphs_tmp.nii -a FMmag_tmp.nii -o FMunw_tmp.nii -m FMmsk_tmp.nii
!gunzip -f FMunw_tmp.nii.gz

tmp = load_nii('FMunw_tmp.nii');
% FM = -1*double(tmp.img)/((rInfo.TE(2)-rInfo.TE(1))*(1e-6));
FM = -1*double(tmp.img)/((rInfo.TE(3)-rInfo.TE(2))*(1e-6));

% if nargout > 1
% %% Expand using smoothn
%     FM(~mask) = NaN;
%     FMsmoothed = smoothn(FM).*maskX;
% end

FM(isnan(FM))=0;

%% Smooth field map over entire image, enforcing data fidelity only over the tight mask used in SENSE Map estimation
beta = 1;
niter = 100;

FMsmoothed = zeros(rInfo.N,rInfo.N,rInfo.nSlices);
TEMask = logical(ones(rInfo.nEchoes,1));

lastTEIndex = find(TEMask,1,'last');
maskSmoothing = maskX;
maskSmoothing = reshape(maskSmoothing,rInfo.N,rInfo.N,rInfo.nSlices);
data_weight = double(maskSmoothing);

R = Robj(logical(ones(size(FMsmoothed(:,:,1)))),'beta',beta,'potential','quad');
A = sensemult(ones(size(FMsmoothed(:,:,1))));

for sliceIndex = 1:rInfo.nSlices
    init_img = FM(:,:,sliceIndex).*mask(:,:,sliceIndex);
    init_data = FM(:,:,sliceIndex).*mask(:,:,sliceIndex);
    W = Gdiag(col(data_weight(:,:,sliceIndex)));
    tmp = solve_pwls_pcg(init_img(:), A,W ,init_data(:), R, 'niter',niter,'stepper',{'qs1',1});
    tmp = reshape(tmp,rInfo.N,rInfo.N);
    FMsmoothed(:,:,sliceIndex) = tmp;
end

FMsmoothed = FMsmoothed.*maskX;
FMsmoothed(isnan(FM))=0;

!rm FMmag_tmp.nii
!rm FMphs_tmp.nii
!rm FMmsk_tmp.nii
!rm FMunw_tmp.nii
