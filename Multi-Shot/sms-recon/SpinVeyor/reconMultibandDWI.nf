#!/usr/bin/env nextflow


/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */
params.workingDir = "$HOME/project"
params.subjectID = "DefaultSubject"
params.senFieldMap = "senFM.dat"
params.MREData = "MRE.dat"

senFMDatFile = file(params.senFieldMap)
mreDatFile = file(params.MREData)
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
    matlab -nodisplay -nodesktop -r "run('$HOME/startup.m'); reconSenFM_Nextflow('!{senFM_Siemens_twix_dat}');"
    """
}


process prepMultibandDWI {

    module "matlab/R2017b"
    cpus 6

    input: 
    file pgse_Siemens_twix_dat from pgseDatFile
    val subjID from subjectID
    file senFM from calibrationData

    output: 
    file "${subjID}.h5" into ISMRMRDFiles

    shell:
    """
    matlab -nodisplay -nodesktop -r "run('$HOME/startup.m'); prepMultibandDWI_Nextflow('!{senFM}','!{pgse_Siemens_twix_dat}','${subjID}');"
    """

}

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
    /opt/PowerGrid/bin/PowerGridPcSense -i !{prepFile} -x 120 -y 120 -z 4 -n 20 -s 4 -D2 -B 1000 -o ./
    """

}

/*
process convert2dNIIsTo4DNIIs {
    
    module : "fsl/5.0.10"

    input : 
    files pcSENSE_Slice*_Rep*_Avg*_Echo*_Phase*_mag.nii into MagNIIs
    files pcSENSE_Slice*_Rep*_Avg*_Echo*_Phase*_phs.nii into PhaseNIIs
    
    output :

}
*/
