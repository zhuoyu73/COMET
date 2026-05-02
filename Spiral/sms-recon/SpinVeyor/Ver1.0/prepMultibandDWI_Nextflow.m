function prepMultibandDWI_Nextflow(fieldMapDataFilename,twixDatFilename, outputFilename)
% function prepMultibandDWI_Nextflow(datadir,filename)
%
% Required Inputs (optional):
%
% - fieldMapDataFilename:  string pointing to MAT file containing field map
%                        data
%
% - twixDatFilename: string with filename with raw data in Siemens TWIX
%                   format
%
% - outputFilename: name for output H5 file in ISMRMRD format
%
% Wrapper function to prep multiband MRE data for reconstruction:
% 2) resamples and reshapes for MB
% 3) recons navigator images
% 4) calculates phase error maps
% 5) converts and saves as ismrmrd format for PowerGrid
%
% Authors:
% Alex Cerjanic - University of Illinois at Urbana-Champaign
% Curtis Johnson - University of Delaware
% Jun 2018
%
% Note: Built out of the former scratch script "reconMultibandMRE" created
% during MB sequence/recon development. Cleaned and consolidated by CLJ.


load(fieldMapDataFilename,'sen','FM','mask','maskC','rInfoSen');

% parse data
rInfo = recoInfo(twixDatFilename);

% Reshape the field map, mask, and sense maps
senResampled    = complex(zeros(rInfo.N,rInfo.N,rInfo.nSlices*rInfo.multibandFactor,rInfo.nCoils));
senResampledNav = complex(zeros(rInfo.NNav,rInfo.NNav,rInfo.nSlices*rInfo.multibandFactor,rInfo.nCoils));

% resample the sense maps for image and nav data
for ii = 1:rInfo.nCoils
    senResampled(:,:,:,ii)    = resampleMap(sen(:,:,:,ii),rInfo,rInfoSen);
    senResampledNav(:,:,:,ii) = resampleMapNav(sen(:,:,:,ii),rInfo,rInfoSen);
end

% resample the field maps for image and nav data
FMResampled       = resampleMap(FM,rInfo,rInfoSen);
FMResampledNav    = resampleMapNav(FM,rInfo,rInfoSen);
%FMImagesResampled = resampleMap(FMImages(:,:,:,2),rInfo,rInfoSen);

% resample the masks for image and nav data
maskResampled    = resampleMap(double(squeeze(mask)),rInfo,rInfoSen);
maskResampled    = maskResampled > 0.5;
maskResampledNav = resampleMapNav(double(squeeze(mask)),rInfo,rInfoSen);
maskResampledNav = maskResampledNav > 0.5;

% multiband pulse phase
% from Wong ISMRM 2012, #2209
lslice_phase = getMbRfPhases(rInfo.multibandFactor);


FMMB         = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices);
%FMImagesMB   = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices);

FMMBNav      = zeros(rInfo.NNav,rInfo.NNav,rInfo.multibandFactor,rInfo.nSlices);

senMB        = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices,rInfo.nCoils);
senMBNav     = zeros(rInfo.NNav,rInfo.NNav,rInfo.multibandFactor,rInfo.nSlices,rInfo.nCoils);

maskMB       = ones(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices);
maskMBNav    = ones(rInfo.NNav,rInfo.NNav,rInfo.multibandFactor,rInfo.nSlices);

% reorder everything for MB format and apply pulse phase
sliceOrder   = 1:rInfo.nSlices*rInfo.multibandFactor;
nn=0;
for jj=1:rInfo.multibandFactor
    for ii=1:rInfo.nSlices
        nn=nn+1;
        FMMB(:,:,jj,ii) = FMResampled(:,:,sliceOrder(nn));
        FMMBNav(:,:,jj,ii) = FMResampledNav(:,:,sliceOrder(nn));
        %FMImagesMB(:,:,jj,ii) = FMImagesResampled(:,:,sliceOrder(nn));
        
        maskMB(:,:,jj,ii) = maskResampled(:,:,sliceOrder(nn));
        maskMBNav(:,:,jj,ii) = maskResampledNav(:,:,sliceOrder(nn));
        for kk=1:rInfo.nCoils
            senMB(:,:,jj,ii,kk) = senResampled(:,:,sliceOrder(nn),kk).*exp(-1i*lslice_phase(jj));
            senMBNav(:,:,jj,ii,kk) = senResampledNav(:,:,sliceOrder(nn),kk).*exp(-1i*lslice_phase(jj));
        end
    end
end

% Flip Slices
senMB = flip(senMB,3);
senMBNav = flip(senMBNav,3);
FMMB = flip(FMMB,3);
FMMBNav = flip(FMMBNav,3);
%FMImagesMB = flip(FMImagesMB,3);

% reconstruct navigator images

imgNav = fieldCorrectedNavRecon(rInfo,senMBNav,maskMBNav,FMMBNav,'Rbeta',1000,'dims2penalize',[1,1,0],'Niter',20,'L',0); % CLJ: turned up the beta and iterations
save imgNav.mat imgNav % To be removed in the future 
% warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
% warning('!!!reconstructing with zero PMaps for debugging!!!');
% warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
% 
% PMaps = zeros(rInfo.N, rInfo.N, rInfo.multibandFactor, rInfo.nShots,rInfo.nPartitions,rInfo.Slices,rInfo.nAverages,rInfo.nPhases,rInfo.nEchoes,rInfo.nRepetitions);

% calculate phase error maps

PMaps = calcDWINavPhaseErrors(rInfo, imgNav); %AMC CHANGE!
save -v7.3 PMaps.mat PMaps

% write data out to ISMRMRD format for PowerGrid
if rInfo.nCoils < 32
convertRecoInfoToIsmrmrd(sprintf('%s.h5',outputFilename),rInfo,permute(senMB,[1 2 3 5 4]), FMMB, PMaps);
else
[rInfoCC,senCC] = compressCoils(rInfo,senMB,'energyLevel',0.95);
convertRecoInfoToIsmrmrd(sprintf('%s.h5',outputFilename),rInfoCC,permute(senCC,[1 2 3 5 4]), FMMB, PMaps);
end

