function convertRecoInfoToIsmrmrdCIVIC(filename,rInfo,sen, FM, PhaseMaps)
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
if nargin > 4
    PMaps = ismrmrd.NDArray(col(PhaseMaps));
    dset.appendArray('PhaseMaps', PMaps);
end

% It is very slow to append one acquisition at a time, so we're going
% to append a block of acquisitions at a time.
% In this case, we'll do it one repetition at a time to show off this
% feature.  Each block has nY aquisitions
acqblock = ismrmrd.Acquisition(rInfo.nShots);

% Set the header elements that don't change
acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = rInfo.shotLength;
acqblock.head.center_sample(:) = 0;
acqblock.head.active_channels(:) = rInfo.nCoils;
acqblock.head.read_dir  = repmat([1 0 0]',[1 rInfo.nShots]);
acqblock.head.phase_dir = repmat([0 1 0]',[1 rInfo.nShots]);
acqblock.head.slice_dir = repmat([0 0 1]',[1 rInfo.nShots]);
acqblock.head.trajectory_dimensions = repmat(4,[1 rInfo.nShots]);


% Loop over the acquisitions, set the header, set the data and append
% For things we aren't dealing with for now.
eco = 1;
phs = 1;
rep = 1;
scanCounter = 0;
for avg = 1:rInfo.nAverages
%for avg = 11
    for slc = 1:rInfo.nSlices
        for par = 1:rInfo.nPartitions
            for shot = 1:rInfo.nShots
                
                % Set the header elements that change from acquisition to the next
                % c-style counting
                acqblock.head.scan_counter(shot) = scanCounter;
                scanCounter = scanCounter + 1;
                acqblock.head.idx.kspace_encode_step_1(shot) = shot-1;
                acqblock.head.idx.kspace_encode_step_2(shot) = par-1;
                acqblock.head.idx.repetition(shot) = 0;
                acqblock.head.idx.slice(shot) = slc - 1;
                
                % Set the flags
                acqblock.head.clearAllFlags(shot);
                
                if shot == 1
                    acqblock.head.setFlag('ACQ_FIRST_IN_ENCODE_STEP1', shot);
                end
                
                if shot == rInfo.nShots
                    acqblock.head.setFlag('ACQ_LAST_IN_ENCODE_STEP1', shot);
                end
                
                if par == 1 && shot == 1
                    acqblock.head.setFlag('ACQ_FIRST_IN_ENCODE_STEP2', shot);
                end
                
                if par == rInfo.nPartitions && shot == rInfo.nShots
                    acqblock.head.setFlag('ACQ_LAST_IN_ENCODE_STEP2', shot);
                end
                
                if slc == 1 && shot == 1 && par == 1
                    acqblock.head.setFlag('ACQ_FIRST_IN_SLICE', shot);
                end
                
                if slc == rInfo.nSlices && par == rInfo.nPartitions && shot == rInfo.nShots
                    acqblock.head.setFlag('ACQ_LAST_IN_SLICE', shot);
                end
                                
                if slc == 1 && shot == 1 && par == 1 && avg == 1
                    acqblock.head.setFlag('ACQ_FIRST_IN_REPETITION', shot);
                end
                
                if slc == rInfo.nSlices && par == rInfo.nPartitions && shot == rInfo.nShots && avg == rInfo.nAverages
                    acqblock.head.setFlag('ACQ_LAST_IN_REPETITION', shot);
                end
                              
                % fill the data
                %acqblock.data{shot} = squeeze(data(:,shot,:,rep));
                acqblock.data{shot} = squeeze(rInfo.dataRead(shot,par,slc,avg,phs,eco,rep));
                % attach the trajectory
                %acqblock.traj{shot} = [kx(nR0*(acqno-1)+1:nR0*acqno),ky(nR0*(acqno-1)+1:nR0*acqno),kz(nR0*(acqno-1)+1:nR0*acqno),t(nR0*(acqno-1)+1:nR0*acqno)].';
                acqblock.traj{shot} = [col(rInfo.kRead(rInfo.dataMask,shot,par,slc,avg,phs,eco,rep)), ...
                    col(rInfo.kPhase(rInfo.dataMask,shot,par,slc,avg,phs,eco,rep)), ...
                    col(rInfo.kSlice(rInfo.dataMask,shot,par,slc,avg,phs,eco,rep)), ...
                    col(rInfo.timingVec(rInfo.dataMask,shot,par,slc,avg,phs,eco,rep))].';
                
                

            end % shot loop
                % Append the acquisition block
                dset.appendAcquisition(acqblock);
        end % partition loop
    end % slice loop
end % avg loop



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
if rInfo.multibandFactor < 1 % 3D SMS case
    header.encoding.encodedSpace.fieldOfView_mm.z = rInfo.sliceThickness*rInfo.multibandFactor;
elseif rInfo.nPartitions < 1 % 3D case
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
header.encoding.encodingLimits.kspace_encoding_step_1.maximum = rInfo.nShots-1;
header.encoding.encodingLimits.kspace_encoding_step_1.center = floor((rInfo.nShots-1)/2);
header.encoding.encodingLimits.kspace_encoding_step_2.minimum = 0;

% Step 1 is partitions (kz encodes)
if rInfo.multibandFactor > 1 % 3D SMS case
   header.encoding.encodingLimits.kspace_encoding_step_2.maximum = rInfo.multibandFactor-1;
   header.encoding.encodingLimits.kspace_encoding_step_2.center = floor((rInfo.multibandFactor-1)/2);
elseif rInfo.nPartitions > 1 % 3D Encoded case
    header.encoding.encodingLimits.kspace_encoding_step_2.maximum = rInfo.nPartitions-1;
    header.encoding.encodingLimits.kspace_encoding_step_2.center = floor((rInfo.nPartitions-1)/2);
else % 2D Case
    header.encoding.encodingLimits.kspace_encoding_step_2.maximum = 0;
    header.encoding.encodingLimits.kspace_encoding_step_2.center = 0;
end

% Deal with Slices
header.encoding.encodingLimits.slice.minimum = 0;
header.encoding.encodingLimits.slice.maximum = rInfo.nSlices-1;
header.encoding.encodingLimits.slice.center = floor((rInfo.nSlices-1)/2);

% Deal with Repetitions
header.encoding.encodingLimits.repetition.minimum = 0;
header.encoding.encodingLimits.repetition.maximum = rInfo.nAverages-1;
%header.encoding.encodingLimits.repetition.maximum = 1;
header.encoding.encodingLimits.repetition.center = floor((rInfo.nAverages-1)/2);
%header.encoding.encodingLimits.repetition.center = 0;

%% Serialize and write to the data set
xmlstring = ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);

%% Write the dataset
dset.close();

end
