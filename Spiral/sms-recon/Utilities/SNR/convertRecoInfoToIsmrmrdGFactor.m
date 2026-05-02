function convertRecoInfoToIsmrmrdGFactor(filename,rInfo,sen, FM, NReps, shotList, parList,sliceList)
%convertRecoInfoToIsmrmrd Adapter from recoInfo object to ismrmrd file
%
%
%
% for SMS MRE: sense maps need to be [ny nx nz ncoil nslice]

% Create an empty ismrmrd dataset
if exist(filename,'file')
    error(['File ' filename ' already exists.  Please remove first'])
    %delete(filename);
end
dset = ismrmrd.Dataset(filename);

% load the coil sensitivities

% create NDArray with the data
SENSEMap = ismrmrd.NDArray(col(sen));

% Write to file
dset.appendArray('SENSEMap',SENSEMap);

% load the field map
FM = col(FM);

% create NDArray with the data
FieldMap = ismrmrd.NDArray(FM);

% Write to file
dset.appendArray('FieldMap',FieldMap);

% if exists append PMaps to file

% It is very slow to append one acquisition at a time, so we're going
% to append a block of acquisitions at a time.
% In this case, we'll do it one repetition at a time to show off this
% feature.  Each block has nY aquisitions
acqblock = ismrmrd.Acquisition(length(shotList));

% Set the header elements that don't change
acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = rInfo.shotLength;
acqblock.head.center_sample(:) = 0;
acqblock.head.active_channels(:) = rInfo.nCoils;
acqblock.head.read_dir  = repmat([1 0 0]',[1 length(shotList)]);
acqblock.head.phase_dir = repmat([0 1 0]',[1 length(shotList)]);
acqblock.head.slice_dir = repmat([0 0 1]',[1 length(shotList)]);
acqblock.head.trajectory_dimensions = repmat(4,[1 length(shotList)]);


% Loop over the acquisitions, set the header, set the data and append
% For things we aren't dealing with for now.
eco = 1;
phs = 1;
avg = 1;
scanCounter = 0;
for rep = 1:NReps
    for slc = 1:length(sliceList)
        for par = 1:1
            for shot = 1:length(shotList)
                
                % Set the header elements that change from acquisition rInfo.nShotsto the next
                % c-style counting
                acqblock.head.scan_counter(shot) = scanCounter;
                scanCounter = scanCounter + 1;
                acqblock.head.idx.kspace_encode_step_1(shot) = shot-1;
                acqblock.head.idx.kspace_encode_step_2(shot) = par-1;
                acqblock.head.idx.repetition(shot) = rep - 1;
                acqblock.head.idx.slice(shot) = slc - 1;
                
                % Set the flags
                acqblock.head.clearAllFlags(shot);
                
                if shot == 1
                    acqblock.head.setFlag('ACQ_FIRST_IN_ENCODE_STEP1', shot);
                end
                
                if shot == length(shotList)
                    acqblock.head.setFlag('ACQ_LAST_IN_ENCODE_STEP1', shot);
                end
                
                if par == 1 && shot == 1
                    acqblock.head.setFlag('ACQ_FIRST_IN_ENCODE_STEP2', shot);
                end
                
                if par == 1 && shot == length(shotList)
                    acqblock.head.setFlag('ACQ_LAST_IN_ENCODE_STEP2', shot);
                end
                
                if slc == 1 && shot == 1 && par == 1
                    acqblock.head.setFlag('ACQ_FIRST_IN_SLICE', shot);
                end
                
                if slc == length(sliceList) && par == 1 && shot == length(shotList)
                    acqblock.head.setFlag('ACQ_LAST_IN_SLICE', shot);
                end
                                
                if slc == 1 && shot == 1 && par == 1 && rep == 1
                    acqblock.head.setFlag('ACQ_FIRST_IN_REPETITION', shot);
                end
                
                
                if slc == length(sliceList) && par == 1 && shot == length(shotList) && rep == NReps
                    acqblock.head.setFlag('ACQ_LAST_IN_REPETITION', shot);
                end
                              
                % fill the data
                %acqblock.data{shot} = squeeze(data(:,shot,:,rep));
                % Genearte the synthetic noise for a Gfactor map calculation via
                % Pseudoreplica method
                temp = squeeze(rInfo.dataRead(shotList(shot),parList(shot),sliceList(slc),avg,phs,eco,1));
                sizesData = size(temp);
                
                if (NReps == 1)
                    noise = 0;
                else
                    noise = mvnrnd(zeros(rInfo.nCoils,1),rInfo.noiseCorr,sizesData(1));
                end
                
                acqblock.data{shot} = temp + noise;
                
                % attach the trajectory
                %acqblock.traj{shot} = [kx(nR0*(acqno-1)+1:nR0*acqno),ky(nR0*(acqno-1)+1:nR0*acqno),kz(nR0*(acqno-1)+1:nR0*acqno),t(nR0*(acqno-1)+1:nR0*acqno)].';
                acqblock.traj{shot} = [col(rInfo.kRead(rInfo.dataMask,shotList(shot),parList(shot),sliceList(slc),avg,phs,eco,1)), ...
                    col(rInfo.kPhase(rInfo.dataMask,shotList(shot),parList(shot),sliceList(slc),avg,phs,eco,1)), ...
                    col(rInfo.kSlice(rInfo.dataMask,shotList(shot),parList(shot),sliceList(slc),avg,phs,eco,1)), ...
                    col(rInfo.timingVec(rInfo.dataMask,shotList(shot),parList(shot),sliceList(slc),avg,phs,eco,1))].';
                
                
            end % shot loop
                % Append the acquisition block
                dset.appendAcquisition(acqblock);
        end % partition loop
    end % slice loop
end % rep loop



%% Fill the xml header
% We create a matlab struct and then serialize it to xml.
% Look at the xml schema to see what the field names should be

header = [];

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = ...
                        str2double(rInfo.DataObject.hdr.Dicom.DICOM.lFrequency); 

% Acquisition System Information (Optional)
header.acquisitionSystemInformation.systemVendor = rInfo.scannerManufacturer;
header.acquisitionSystemInformation.systemModel = rInfo.scannerModel;
header.acquisitionSystemInformation.receiverChannels = rInfo.nCoils;

% The Encoding (Required)
header.encoding.trajectory = 'spiral'; % Probably a safe bet for our lab
header.encoding.encodedSpace.fieldOfView_mm.x = rInfo.FOV;
header.encoding.encodedSpace.fieldOfView_mm.y = rInfo.FOV;
if rInfo.multibandFactor > 1 % 3D SMS case
    header.encoding.encodedSpace.fieldOfView_mm.z = rInfo.sliceThickness*rInfo.multibandFactor;
elseif rInfo.nPartitions > 1 % 3D case
    header.encoding.encodedSpace.fieldOfView_mm.z = rInfo.sliceThickness*rInfo.nPartitions;
else % 2D case
    header.encoding.encodedSpace.fieldOfView_mm.z = rInfo.sliceThickness;
end

header.encoding.encodedSpace.matrixSize.x = rInfo.N;
header.encoding.encodedSpace.matrixSize.y = rInfo.N;

if rInfo.multibandFactor > 1
    header.encoding.encodedSpace.matrixSize.z = rInfo.multibandFactor;
elseif rInfo.nPartitions > 1
    header.encoding.encodedSpace.matrixSize.z = rInfo.nPartitions;
else
    header.encoding.encodedSpace.matrixSize.z = 1;
end

% Recon Space
% (in this case same as encoding space)
header.encoding.reconSpace = header.encoding.encodedSpace;
% Encoding Limits;

% Step 0 is shots for us
header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_1.maximum = length(shotList)-1;
header.encoding.encodingLimits.kspace_encoding_step_1.center = floor((length(shotList)-1)/2);
header.encoding.encodingLimits.kspace_encoding_step_2.minimum = 0;

% Step 1 is partitions (kz encodes)
if rInfo.nPartitions > 1 % 3D Encoded case
    header.encoding.encodingLimits.kspace_encoding_step_2.maximum = 0;
    header.encoding.encodingLimits.kspace_encoding_step_2.center = 0;
else % 2D Case
    header.encoding.encodingLimits.kspace_encoding_step_2.maximum = 0;
    header.encoding.encodingLimits.kspace_encoding_step_2.center = 0;
end

% Deal with Slices
header.encoding.encodingLimits.slice.minimum = 0;
header.encoding.encodingLimits.slice.maximum = length(sliceList)-1;
header.encoding.encodingLimits.slice.center = floor((length(sliceList)-1)/2);

% Deal with Repetitions
header.encoding.encodingLimits.repetition.minimum = 0;
header.encoding.encodingLimits.repetition.maximum = NReps-1;
header.encoding.encodingLimits.repetition.center = floor((NReps-1)/2);

%% Serialize and write to the data set
xmlstring = ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);

%% Write the dataset
dset.close();

end

