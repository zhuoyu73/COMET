# Pipelines for Iterative Reconstruction

This folder contains a number of [Nextflow](http://nextflow.io) pipelines for reconstructing images in an automated fashion. Nextflow separates the differences between implementation of the pipeline steps from the resource management of various clusters, clouds, or local machines. Ideally, this allows the pipleline maintainer to focus on the implementation of each step in an abstracted (and portable) fashion.

## Piplelines
- reconMultibandDWI.nf - Multiband DWI recon as appropriate for DTI or IVIM.
