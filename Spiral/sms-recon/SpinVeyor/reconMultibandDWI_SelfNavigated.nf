#!/usr/bin/env nextflow


/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */
params.workingDir = "$HOME/project"
params.subjectID = "DefaultSubject"
params.senFieldMap = "senFM.dat"
params.DWIData = "DWI.dat"

senFMDatFile = file(params.senFieldMap)
dwiDatFile = file(params.DWIData)
subjectID = Channel.value(params.subjectID)

//    file PGSE_Siemens_twix_dat from datFiles


process reconSenFM {
    
    module "matlab/R2017b"
    cpus 6

    input:
    file senFM_Siemens_twix_dat from senFMDatFile
    
    output: 
    file 'senFM.mat' into calibrationData

    shell:
    """
    matlab -nodisplay -nodesktop -r "run('/shared/mrfil-data/acerja2/repos/test/startup.m'); reconSenFM_Nextflow('!{senFM_Siemens_twix_dat}');"
    """
}


process prepMultibandDWI_SelfNavigated {

    module "matlab/R2017b"
    cpus 6

    input: 
    file dwi_Siemens_twix_dat from dwiDatFile
    val subjID from subjectID
    file senFM from calibrationData

    output: 
    file 'Shot*.h5' into ISMRMRDFiles

    shell:
    """
    matlab -nodisplay -nodesktop -r "run('/shared/mrfil-data/acerja2/repos/test/startup.m'); prepMultibandDWI_SelfNavigated_Nextflow('!{senFM}','!{dwi_Siemens_twix_dat}','${subjID}');"
    """

}

process runPGReconNavigators {

    container = "powergrid"
    containerOptions = "--runtime=nvidia "
	maxFork = 1

    input:
    val subjID from subjectID
    each file(prepFile) from ISMRMRDFiles
    
    output: 
    file 'pcSENSE_Slice*_Rep*_Avg*_Echo*_Phase*_mag.nii' into MagNIIs
    file 'pcSENSE_Slice*_Rep*_Avg*_Echo*_Phase*_phs.nii' into PhaseNIIs

    shell:
    """
	export OMP_NUM_THREADS=2;
    /opt/PowerGrid/bin/PowerGridIsmrmrd -i !{prepFile} -I hanning -t 20 -F DFT -x 120 -y 120 -z 2 -n 30 -s 1 -D 2 -B100 -o ./
    """

}

/*
process runPGRecon {
    container = "acerja2/privatepg"
    containerOptions = "--runtime=nvidia "


    input:
    val subjID from subjectID
    file prepFile from ISMRMRDFiles
    
    output: 
    file 'pcSENSE_Slice*_Rep*_Avg*_Echo*_Phase*_mag.nii' into MagNIIs
    file 'pcSENSE_Slice*_Rep*_Avg*_Echo*_Phase*_phs.nii' into PhaseNIIs

    shell:
    """
    /opt/PowerGrid/bin/PowerGridPcSense -i !{prepFile} -x 120 -y 120 -z 4 -n 20 -s 2 -D2 -B 1000 -o ./
    """

}


process convert2dNIIsTo4DNIIs {
    
    module : "fsl/5.0.10"

    input : 
    files pcSENSE_Slice*_Rep*_Avg*_Echo*_Phase*_mag.nii into MagNIIs
    files pcSENSE_Slice*_Rep*_Avg*_Echo*_Phase*_phs.nii into PhaseNIIs
    
    output :

}
*/
