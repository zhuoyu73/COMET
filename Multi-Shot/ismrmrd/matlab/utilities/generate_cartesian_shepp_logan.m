function generate_cartesian_shepp_logan(filename, varargin)


%% Parse Inputs
p = inputParser;
defaultMatrix = 256;
defaultCoils = 8;
defaultReadoutOversampling = 2;
defaultRepetitions = 1;
defaultAccelerationFactor = 1;
defaultNoiseLevel = 0.05;
defaultDatasetName = 'dataset';
defaultSaveNoiseCalibration = false;
defaultSaveKSpace = false;

addOptional(p,'matrixSize',defaultMatrix,@isnumeric);
addOptional(p,'numCoils',defaultCoils,@isnumeric);
addOptional(p,'readoutOversampling',defaultReadoutOversampling, @isnumeric);
addOptional(p,'accelerationFactor',defaultAccelerationFactor,@isnumeric);
addOptional(p,'noiseLevel',defaultNoiseLevel,@isnumeric);
addOptional(p,'numRepetitions',defaultRepetitions,@isnumeric);
addOptional(p,'saveKSpace',defaultSaveKSpace,@isboolean);
addOptional(p,'saveNoiseCalibration',defaultSaveNoiseCalibration,@isboolean);
addOptional(p,'datasetName',defaultDatasetName,@isstring);

parse(p,varargin{:});

% Prepare variables for use
matrixSize = p.Results.matrixSize;
nCoils = p.Results.numCoils;
ros = p.Results.readoutOversampling;
RFactor = p.Results.accelerationFactor;
noiseLevel = p.Results.noiseLevel;
nReps = p.Results.numRepetitions;
saveKSpace = p.Results.saveKSpace;
saveNoiseCal = p.Results.saveNoiseCalibration;
datasetName = p.Results.datasetName;

%% Generate phantom object and coil sensitivies

% Use built in Matlab function
phantomObject = phantom('Modified Shepp-Logan',matrixSize);

% Generate coil sensitivies
csm = generate_birdcage_sensitivities(matrixSize,nCoils,1.5);

% Generate coil images
coilImages = zeros(matrixSize*matrixSize,nCoils);
for cc = 1:nCoils
    coilImages(:,cc) = col(csm(:,:,cc)).*col(phantomObject);
end
coilImages = reshape(coilImages,[matrixSize,matrixSize,nCoils]);

%% Begin creating ismrmrd dataset (and file)

% Create an empty ismrmrd dataset
if exist(filename,'file')
    error(['File ' filename ' already exists.  Please remove first'])
end

dset = ismrmrd.Dataset(filename);

readoutLength = ros*matrixSize; % Account for readout oversampling


if saveNoiseCal
    acq = ismrmrd.Acquisition();
    acq.head.setFlag('ACQ_IS_NOISE_MEASUREMENT',1);
    acq.head.sample_time_us(:) = 5.0;
    %Generate noise for noise calibration data
    noiseData = randn(readoutLength,ncoils) + 1j*randn(readoutLength,ncoils);
    
    acq.data{1} = noiseData;
    dset.appendAcquisition(acq);
end

if saveKSpace
    %acq = ismrmrd.Acquisition(readoutLength,nCoils,2);
    acq = ismrmrd.Acquisition(); % Don't need to prespecify size in Matlab
else
    %acq = ismrmrd.Acquisition(readoutLength,nCoils);
    acq = ismrmrd.Acquisition(); % Don't need to prespecify size in Matlab
end

acq.head.available_channels(:) = nCoils;
acq.head.center_sample(:) = floor(readoutLength/2);
acq.head.number_of_samples(:) = readoutLength;

% Synthesize the fully sampled k-space data
fullKSpaceRep = zeros(readoutLength, matrixSize, nCoils, nReps);
for rep = 1:nReps
    for coil = 1:nCoils
        noise = noiseLevel * (randn(readoutLength,matrixSize)+1j*randn(readoutLength,matrixSize));
        tempCoilImage = zeros(readoutLength,matrixSize); %Deal with readout oversampling
        tempCoilImage(ceil(readoutLength/4):(3*floor(readoutLength/4)-1),:) = coilImages(:,:,coil);
        fullKSpaceRep(:,:,coil,rep) = 1/numel(tempCoilImage)*fftshift(fft2(fftshift(tempCoilImage + noise)));
    end
end
for kk = 1:nReps
    % Deal with cartesian undersampling
    for ii = 1:RFactor
        
        for jj = ii:RFactor:matrixSize % Phase encode loop
            
            acq.head.clearAllFlags();
            if (jj == ii)
                acq.head.setFlag('ACQ_FIRST_IN_SLICE');
            end
            
            if (ii >= (matrixSize-RFactor))
                acq.head.setFlag('ACQ_LAST_IN_SLICE');
            end
            acq.head.idx.kspace_encode_step_1 = jj-1;
            acq.head.idx.repetition = kk-1;
            acq.head.sample_time_us = 5.0;
            acq.data{1} = single(squeeze(fullKSpaceRep(:,jj,:,kk)));
            acq.head.active_channels = nCoils;
            
            if saveKSpace
                ky = ii - floor(matrixSize/2)/matrixSize;
                kx = -floor(matrixSize/2):1:floor(matrixSize/2);
                acq.traj = horzcat(kx(:),ky(:));
            end
            
            dset.appendAcquisition(acq);
        end
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill the xml header %
%%%%%%%%%%%%%%%%%%%%%%%%
% We create a matlab struct and then serialize it to xml.
% Look at the xml schema to see what the field names should be

header = [];

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = 63500000; % 3T

% Acquisition System Information (Optional)
header.acquisitionSystemInformation.systemVendor = 'ISMRM Synthetic Imaging Lab';
header.acquisitionSystemInformation.receiverChannels = nCoils;

% The Encoding (Required)
header.encoding.trajectory = 'cartesian';
header.encoding.encodedSpace.fieldOfView_mm.x = 600;
header.encoding.encodedSpace.fieldOfView_mm.y = 300;
header.encoding.encodedSpace.fieldOfView_mm.z = 6;
header.encoding.encodedSpace.matrixSize.x = readoutLength;
header.encoding.encodedSpace.matrixSize.y = matrixSize;
header.encoding.encodedSpace.matrixSize.z = 1;

% Recon Space
header.encoding.reconSpace.matrixSize.x = matrixSize;
header.encoding.reconSpace.matrixSize.y = matrixSize;
header.encoding.reconSpace.matrixSize.z = 1;
header.encoding.reconSpace.fieldOfView_mm.x = 300;
header.encoding.reconSpace.fieldOfView_mm.y = 300;
header.encoding.reconSpace.fieldOfView_mm.z = 6;

% Encoding Limits
% Replacement for Limit struct
% Limit(minimum, maximum, center) => Limit.minimum; Limit.maximum; Limit.center;

header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_0.maximum = matrixSize-1;
header.encoding.encodingLimits.kspace_encoding_step_0.center = floor(matrixSize/2);
% Replacement for Limit struct
% Limit(minimum, maximum, center) => Limit.minimum; Limit.maximum; Limit.center;
header.encoding.encodingLimits.repetition.minimum = 0;
header.encoding.encodingLimits.repetition.maximum = nReps*RFactor-1;
header.encoding.encodingLimits.repetition.center = 0;

if (RFactor > 1)
    header.encoding.parallel.accelerationFactor.kspace_encoding_step_1 = RFactor;
    header.encoding.parallel.accelerationFactor.kspace_encoding_step_2 = 1;
    header.encoding.parallel.calibrationMode = 'interleaved';
end

%% Serialize and write to the data set
xmlstring = ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);

% Save the phantom as an NDArray
phantomNDArray = ismrmrd.NDArray(phantomObject);
dset.appendArray('phantom',phantomNDArray);

% Save the csm as an NDArray
csmNDArray = ismrmrd.NDArray(csm);
dset.appendArray('csm',csmNDArray);

% Save the coil images as an NDArray
coilImagesNDArray = ismrmrd.NDArray(coilImages);
dset.appendArray('coilImages',coilImagesNDArray);

%% Close the dataset
dset.close();

end

function csm = generate_birdcage_sensitivities(matrixSize, nCoils, relativeRadius)

csm = zeros(matrixSize, matrixSize, nCoils);
for kk = 1:matrixSize % x
    for jj = 1:matrixSize % y
        for ii = 1:nCoils % coil index
            coilX = relativeRadius*cos(ii*(2*pi/nCoils));
            coilY = relativeRadius*cos(ii*(2*pi/nCoils));
            coilPhase = -ii*(2*pi/nCoils);
            yCo = (jj - 2*matrixSize)/(matrixSize/2) - coilY;
            xCo = (kk - 2*matrixSize)/(matrixSize/2) - coilX;
            rr = sqrt(xCo.^2 + yCo.^2);
            phi = atan2(xCo,-yCo) + coilPhase;
            csm(kk,jj,ii) = 1./rr.*exp(1j*phi);
        end
    end
end
end