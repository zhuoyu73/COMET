% Wrapper script to prepare CIVIC reconstructions for PowerGrid
% Run from the directory with the PGSE or STEAM data.
curFolder = pwd;

cd senfm
rInfoSen = recoInfo();

if ~exist('FM.mat','file')
   reconSenFMCIVIC;
end

load FM.mat
load FMImages.mat
load sen.mat
load mask.mat

cd ..

% parse data
rInfo = recoInfo();

% Reshape the field map, mask, and sense maps
senResampled    = complex(zeros(rInfo.N,rInfo.N,rInfo.nPartitions,rInfo.nCoils));
senResampledNav = complex(zeros(rInfo.NNav,rInfo.NNav,rInfo.nPartitionsNav,rInfo.nCoils));

for ii = 1:rInfo.nCoils
   % senResampled(:,:,:,ii)    = resample_map_resolution(sen(:,:,:,ii),rInfo.N,rInfo.nSlices*rInfo.multibandFactor,rInfoSen.FOV*10,rInfo.FOV*10);
   % senResampledNav(:,:,:,ii) = resample_map_resolution(sen(:,:,:,ii),rInfo.NNav,rInfo.nSlices*rInfo.multibandFactor,rInfoSen.FOV*10,rInfo.FOV*10);
   senResampled(:,:,:,ii) =    resampleMap(sen(:,:,:,ii),rInfo,rInfoSen);
   senResampledNav(:,:,:,ii) = resampleMapNav(sen(:,:,:,ii),rInfo,rInfoSen);
end

%FMResampled      = resample_map_resolution(FM,rInfo.N,rInfo.nSlices*rInfo.multibandFactor,rInfoSen.FOV*10,rInfo.FOV*10);
%FMResampledNav   = resample_map_resolution(FM,rInfo.NNav,rInfo.nSlices*rInfo.multibandFactor,rInfoSen.FOV*10,rInfo.FOV*10);
FMResampled      = resampleMap(FM,rInfo,rInfoSen);
FMResampledNav   = resampleMapNav(FM,rInfo,rInfoSen);
%FMImagesResampled      = resampleMap(FMImages(:,:,:,2),rInfo,rInfoSen);

maskResampled    = resampleMap(double(squeeze(mask)),rInfo,rInfoSen);
maskResampled = maskResampled > 0.5;
%maskResampleNav = true(rInfo.NNav,rInfo.NNav,rInfo.nPartitionsNav);
maskResampledNav = resampleMapNav(double(squeeze(mask)),rInfo,rInfoSen);
maskResampledNav = maskResampledNav > 0.5;

senResampledNav = permute(senResampledNav,[1,2,3,5,4]);
%FMResampledNav = permute(FMResampledNav,[1,2,3,4]);
maskResampledNav = true(40,40,20);

%senResampledNav = flip(senResampledNav,2);
%senResampledNav = flip(senResampledNav,3);
%FMResampledNav = flip(FMResampledNav,3);
rInfo.dataMask = true(rInfo.shotLength,1);
cImages = gridCoilImages(rInfo);
save -v7.3 cImages.mat cImages
% Calculate DCF for navigator
% kRead = col(rInfo.kReadNav(rInfo.dataMaskNav,1,1,1,1,1,1,1));
% kPhase = col(rInfo.kPhaseNav(rInfo.dataMaskNav,1,1,1,1,1,1,1));
% kSlice = col(rInfo.kSliceNav(rInfo.dataMaskNav,1,1,1,1,1,1,1));
% 
% dcf = calc3DDensityCompenasation(kRead,kPhase,kSlice);

% if ~exist('imgNav.mat','file')
%     imgNav = fieldCorrectedNavReconCIVIC(rInfo,senResampledNav,maskResampledNav,FMResampledNav, dcf,'L',0,'Rbeta',200,'dims2penalize',[1,1,1],'Niter',10);
%     save imgNav.mat imgNav
% else
%     load imgNav.mat
% end

%for mm = 1:rInfo.nAverages
% mm = 6;
%     PMaps = zeros(rInfo.N,rInfo.N,rInfo.nPartitions,rInfo.nShots,rInfo.nPartitions,rInfo.nSlices,rInfo.nPhases,rInfo.nEchoes,rInfo.nRepetitions);
%     for pp = 1:rInfo.nRepetitions
%         for oo = 1:rInfo.nEchoes
%             for nn = 1:rInfo.nPhases
%                 for kk = 1:rInfo.nSlices
%                     for jj = 1:rInfo.nPartitions
%                         for ii = 1:rInfo.nShots
%                         	PMaps(:,:,:,ii,jj,kk,nn,oo,pp) = angle(squeeze(resampleMapNavUp(imgNav(:,:,:,ii,jj,kk,mm,nn,oo,pp),rInfo,rInfo)));
%                         end
%                     end
%                 end
%             end
%         end
%     end
% %end
                
%img = fieldCorrectedReconMB(rInfo, senMB, maskMB, FMMB,'Rbeta',1,'dims2penalize',[1,1,0]);
%senResampled = flip(senResampled,3);

%senResampled = permute(senResampled,[1,2,3,5,4]);

%convertRecoInfoToIsmrmrdCIVIC('data.h5',rInfo,senResampled,FMResampled,PMaps);

%senResampled = flip(senResampled,3);
%maskResampled = flip(maskResampled,3);
%FMResampled = flip(FMResampled,3);
%[rInfoCC, senResampledCC] = compressCoils(rInfo,senResampled,'coilRank',4);
%img = phaseCorrectedReconCIVIC(rInfoCC, senResampledCC, maskResampled, FMResampled, imgNav,'L',0,'averagesToRecon',6,'phasesToRecon',1,'Niter',3);
%img = fieldCorrectedReconMB(rInfo, senMB, maskMB, FMMB,'Rbeta',1,'dims2penalize',[1,1,0]);
%img = phaseCorrectedRecon(rInfo, senMB, maskMB, FMMB, PMaps,'Niter',10,'L',0,'Rbeta',1,'dims2penalize',[1,1,0]);
