function recon_cartesian_2d(filename)

% Create an empty ismrmrd dataset
if isempty(filename)
    error('Specify input filename for ismrmrd H5 file');
end

%% Open the dataset from the file

dset = ismrmrd.Dataset(filename);
hdr = ismrmrd.xml.deserialize(dset.readxml);

% if isempty(hdr.version())
%     disp('XML Header unspecified version.');
% else
%     disp(['XML Header version: ' num2str(hdr.version())]);
% end

if (length(hdr.encoding) ~= 1)
    error('This simple recon only supports one encoding space');
end

nAcq = dset.getNumberOfAcquisitions();

if hdr.encoding(1).encodedSpace.matrixSize.z ~= 1
    error('This simple reconstruction application only supports 2D encoding spaces');
end

nX = hdr.encoding(1).encodedSpace.matrixSize.x;
nY = hdr.encoding(1).encodedSpace.matrixSize.y;
nZ = hdr.encoding(1).encodedSpace.matrixSize.z;

reconNX = hdr.encoding(1).reconSpace.matrixSize.x;
reconNY = hdr.encoding(1).reconSpace.matrixSize.y;
reconNZ = hdr.encoding(1).reconSpace.matrixSize.z;
% To get number of coils, read the first acquisition of data and assume
% that number holds
acq = dset.readAcquisition(1);
nCoils = acq.head.active_channels();

for ii = 1:nAcq
    acq = dset.readAcquisition(ii);
    data(:,ii,:) = acq.data{1};
    
end

disp(['Encoding Matrix Size        : [' num2str(nX) ', ' num2str(nY) ', ' num2str(nZ) ']']);
disp(['Reconstruction Matrix Size  : [' num2str(reconNX) ', ' num2str(reconNY) ', ' num2str(reconNZ) ']'])
disp(['Number of Channels          : ' num2str(nCoils)']);
disp(['Number of Acquisitions      : ' num2str(nAcq) ]);

for ii = 1:nCoils
    coilImagesOS(:,:,ii) = fftshift(fft2(fftshift(data(:,:,ii))));
end

% Remove readout oversampling
offset = floor((nX - reconNX)/2);
coilImages = coilImagesOS(offset+1:offset+reconNX,:,:);

% Take square root of sum of squares image
imageOut = sqrt(sum(abs(coilImages).^2,3));

% Save the image as an ISMRMRD Image
% Klugy workaround to ensure that object is empty 
imageOutObject = ismrmrd.Image();
delete(imageOutObject);
imageOutObject = ismrmrd.Image();


% Add extra guidance as in C++ example
imageOutObject.head.setImageType('MAGNITUDE',1);
imageOutObject.head.slice = 1;
imageOutObject.head.field_of_view = [hdr.encoding(1).reconSpace.fieldOfView_mm.x; ...
    hdr.encoding(1).reconSpace.fieldOfView_mm.y; ...
    hdr.encoding(1).reconSpace.fieldOfView_mm.z];
imageOutObject.data{1} = imageOut;
dset.appendImage('matlab',imageOutObject);

%% Close the dataset
dset.close();

end