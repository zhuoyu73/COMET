function reconSenFM_Nextflow(twixFileInput)
%reconSenFM - Creates SENSE Maps and Field Maps with default options
%
% Syntax:  [cImages, sen, FM, FMImages] = ReconSenFM()
%
% Inputs:
%    None
%
% Outputs:
%    cImages - Gridded Coil Images of "standard dimensions"
%    sen     - Sense Maps
%    FM      - Field Maps
%    FMImages- SENSE Reconstructed images used to compute the field map
%
% Example:
%    [cImages, sen, FM, FMImages] = ReconSenFM();
%    rInfo  = rInfo(filename);
%    images = gridCoilImages(rInfo);
%    im(sum(abs(images).^2,10)); % Sum of Squares image for
%                                            % first echo time in field map
%
%
% Other m-files required: Too numerous to listNew FolderNew Folder
% Subfunctions: none
% MAT-files required: none
%
% Author: Alex Cerjanic
% University of Illinois at Urbana-Champaign
% email address:
% Website:
% 11-Apr-2017; Last revision: 12-Sep-2017


%% Deal with optional input arguments
% Nothing to add at the current time.
% p = inputParser();
TEMask = [];
UseGIRF = false;

% addOptional(p,'UseGIRF',UseGIRF);
% addOptional(p,'TEMask',TEMask);
% parse(p,varargin{:});
% UseGIRF = p.Results.UseGIRF;
% TEMask = p.Results.TEMask;



%% Moving on to actual reconstruction and estimation.
rInfo = recoInfo(twixFileInput,'UseGIRF',UseGIRF);

if isempty(TEMask)
    TEMask = logical(ones(1,rInfo.nEchoes));
end

cImages = gridCoilImages(rInfo);


    if rInfo.nPartitions > 1
        %3D/Multislab style acquisition
        sen = zeros(rInfo.N, rInfo.N, rInfo.nPartitions, rInfo.nSlices,rInfo.nCoils);
        mask = zeros(rInfo.N, rInfo.N, rInfo.nPartitions, rInfo.nSlices,rInfo.nCoils);
        
        for ii = 1:rInfo.nSlices
            senseStack = reshape(cImages(:,:,:,ii,1,1,1,1,1,:),[rInfo.N,rInfo.N,rInfo.nPartitions,rInfo.nCoils]);
            [sen(:,:,:,ii,:), mask(:,:,:,ii,:)] = createSenMap(senseStack, 2);
        end
    else
        %2D Style Acquisition
        senseStack = reshape(cImages(:,:,1,:,1,1,1,1,1,:),[rInfo.N,rInfo.N,rInfo.nSlices,rInfo.nCoils]);
        [sen, mask] = createSenMap(senseStack, 2);
    end

[x, y] = meshgrid(1:rInfo.N,1:rInfo.N);
maskCirc = sqrt((x-(rInfo.N+1)/2).^2+(y-(rInfo.N+1)/2).^2)<((rInfo.N)/2);
maskCirc = repmat(maskCirc,[1 1 rInfo.nPartitions rInfo.nSlices]);
maskC = squeeze(maskCirc);

if rInfo.nPartitions > 1
    % 3D/Multislab style acquisition - To Do: make sure that createFieldMap
    % will process one slab at a time, pass index
    FM = zeros(rInfo.N, rInfo.N, rInfo.nPartitions, rInfo.nSlices);
    FMImages =  zeros(rInfo.N, rInfo.N, rInfo.nPartitions, rInfo.nSlices, rInfo.nEchoes);

    for ii = 1:rInfo.nSlices
        [FM(:,:,:,ii), FMImages(:,:,:,ii,:)] = createFieldMap(rInfo, sen(:,:,:,ii,:), maskCirc(:,:,:,ii,:), 1);
    end
    

        
else
    % FM = zeros(rInfo.N, rInfo.N, rInfo.nSlices);
    %else FMImages = zeros(rInfo.N, rInfo.N, rInfo.nSlices, rInfo.nEchoes);
    [FM, FMImages,FMmask, FMsmoothed] = createFieldMap(rInfo, sen, squeeze(maskCirc(:,:,1,:)), 2, 'nIterations',1,'TEMask',TEMask);
    
    if rInfo.nEchoes == 2 % currently a workaround to allow both Funai and PRELUDE versions
        Fmask = FM~=0;
        [~,FMx] = estimFieldMap_alt(rInfo,FMImages,Fmask);
        FM = FMx;
    end
end 
rInfoSen = rInfo;
save senFM.mat -v7.3 sen mask maskC FM rInfoSen

end

