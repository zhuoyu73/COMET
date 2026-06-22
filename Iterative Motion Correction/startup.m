function startup

setenv( 'FSLDIR', '/usr/local/fsl');
setenv( 'FSLOUTPUTTYPE', 'NIFTI_GZ');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;

addpath(sprintf('%s/sms-recon',pwd));
initializePaths
addpath(sprintf('%s/sms-recon/McIlvain_Routines_seq',pwd));

addpath(genpath(sprintf('%s/ismrmrd/matlab',pwd)))
