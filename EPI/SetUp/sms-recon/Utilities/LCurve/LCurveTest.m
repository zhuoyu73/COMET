% LCurve Test Script - Alex Cerjanic 9/21/2017

%% Establish the reconstruction 

rInfo = recoInfo();

cd senfm
rInfoSen = recoInfo();
if ~exist('FM.mat','file')
   reconSenFM;
end
load FM.mat
load FMImages.mat
load sen.mat
load mask.mat

cd ..

% parse data
rInfo = recoInfo();

% Reshape the field map, mask, and sense maps
senResampled    = complex(zeros(rInfo.N,rInfo.N,rInfo.nSlices*rInfo.multibandFactor,rInfo.nCoils));
senNavResampled = complex(zeros(rInfo.NNav,rInfo.NNav,rInfo.nSlices*rInfo.multibandFactor,rInfo.nCoils));

for ii = 1:rInfo.nCoils
   % senResampled(:,:,:,ii)    = resample_map_resolution(sen(:,:,:,ii),rInfo.N,rInfo.nSlices*rInfo.multibandFactor,rInfoSen.FOV*10,rInfo.FOV*10);
   % senResampledNav(:,:,:,ii) = resample_map_resolution(sen(:,:,:,ii),rInfo.NNav,rInfo.nSlices*rInfo.multibandFactor,rInfoSen.FOV*10,rInfo.FOV*10);
   senResampled(:,:,:,ii) =    resampleMap(sen(:,:,:,ii),rInfo,rInfoSen);
   %senResampledNav(:,:,:,ii) = resampleMap(sen(:,:,:,ii),rInfo,rInfoSen);
end

%FMResampled      = resample_map_resolution(FM,rInfo.N,rInfo.nSlices*rInfo.multibandFactor,rInfoSen.FOV*10,rInfo.FOV*10);
%FMResampledNav   = resample_map_resolution(FM,rInfo.NNav,rInfo.nSlices*rInfo.multibandFactor,rInfoSen.FOV*10,rInfo.FOV*10);
%FMImagesResampled      = resample_map_resolution(FMImages(:,:,:,1),rInfo.N,rInfo.nSlices*rInfo.multibandFactor,rInfoSen.FOV*10,rInfo.FOV*10);
FMResampled      = resampleMap(FM,rInfo,rInfoSen);
%FMResampledNav   = resample_map_resolution(FM,rInfo.NNav,rInfo.nSlices*rInfo.multibandFactor,rInfoSen.FOV*10,rInfo.FOV*10);
FMImagesResampled      = resampleMap(FMImages(:,:,:,1),rInfo,rInfoSen);

maskResampled    = resample_map_resolution(squeeze(mask),rInfo.N,rInfo.nSlices*rInfo.multibandFactor,rInfoSen.FOV*10,rInfo.FOV*10);
%maskResampledNav = resample_map_resolution(squeeze(mask),rInfo.NNav,rInfo.nSlices*rInfo.multibandFactor,rInfoSen.FOV*10,rInfo.FOV*10);

% this is for multiband!! 
lslice_phase = getMbRfPhases(rInfo.multibandFactor);

% reshape sen and FM and mask
FMMB         = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices);
FMImagesMB   = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices);

FMMBNav      = zeros(rInfo.NNav,rInfo.NNav,rInfo.multibandFactor,rInfo.nSlices);

senMB        = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices,rInfo.nCoils);
senMBNav     = zeros(rInfo.NNav,rInfo.NNav,rInfo.multibandFactor,rInfo.nSlices,rInfo.nCoils);

maskMB       = ones(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices);
maskMBNav    = ones(rInfo.NNav,rInfo.NNav,rInfo.multibandFactor,rInfo.nSlices);
%sliceOrder  = [1:2:rInfo.nSlices*rInfo.multibandFactor,2:2:rInfo.nSlices*rInfo.multibandFactor];
sliceOrder   = 1:rInfo.nSlices*rInfo.multibandFactor;
nn=0;
for jj=1:rInfo.multibandFactor
    for ii=1:rInfo.nSlices
        nn=nn+1;
        FMMB(:,:,jj,ii) = FMResampled(:,:,sliceOrder(nn));
        %FMMBNav(:,:,jj,ii) = FMResampledNav(:,:,sliceOrder(nn));
        FMImagesMB(:,:,jj,ii) = FMImagesResampled(:,:,sliceOrder(nn));

        maskMB(:,:,jj,ii) = maskResampled(:,:,sliceOrder(nn));
        for kk=1:rInfo.nCoils
            senMB(:,:,jj,ii,kk) = senResampled(:,:,sliceOrder(nn),kk).*exp(-1i*lslice_phase(jj));
            %senMBNav(:,:,jj,ii,kk) = senResampledNav(:,:,sliceOrder(nn),kk).*exp(-1i*lslice_phase(jj));
            %senMB(:,:,jj,ii,kk) = senResampled(:,:,nn,kk);
        end
    end
end

%senMB       = flip(senMB,3);
%maskMB      = flip(maskMB,3);
%FMMB        = flip(FMMB,3);
%FMImagesMB  = flip(FMImagesMB,3);

%senNav      = flip(senMBNav,3);
%maskNav     = flip(maskMBNav,3);
%FMNav       = flip(FMMBNav,3);

save sen.mat  senMB  senMBNav
save FM.mat   FMMB   FMMBNav
save mask.mat maskMB maskMBNav

rInfo.dataMask = true(rInfo.ShotLength,1);

reconFunction = @(lambda) fieldCorrectedReconMB(rInfo, senMB, maskMB, FMMB, 'Niter', 30, 'L', 0, 'Rbeta', lambda, 'slicesToRecon', 4, 'repetitionsToRecon', 4); 

decadesLambda = -5:1:5;

[ resid, roughness, images ] = calcLCurve(reconFunction, 10.^decadesLambda);

plotLCurve(resid,roughness,10.^decadesLambda);