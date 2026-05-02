# COMET Motion and Wave Field Correction

This repository provides MATLAB pipelines for **COMET motion correction** and **wave field correction** for brain Magnetic Resonance Elastography (MRE) data. The workflow supports both **spiral multiband MRE** and **EPI MRE** acquisitions.

The goal of the pipeline is to:

1. Estimate motion from navigator or image-domain data.
2. Incorporate the estimated motion into image reconstruction or post-processing.
3. Rotate displacement vector fields into a consistent anatomical frame.
4. Generate motion-corrected and wave-field-corrected MRE displacement maps ready for MRE inversion.

---

## Repository Structure

```text
COMET-motion-wavefield-correction/
├── Spiral/
│   ├── startup.m
│   └── sms-recon/
│       ├── initializePaths.m
│       ├── Routines/
│       │   ├── prepMultibandMRE.m
│       │   ├── calc2Dregistrations_ZS.m
│       │   └── reconMultibandMRE_ZS.m
│       └── Utilities/
│           └── MRE/
│               ├── proc_mbmre_1_ZS_RotateVectors.m
│               ├── manual_masking.m
│               └── proc_mbmre_2_ZS_RotateVectors.m
│
└── EPI/
    └── Siemens_EPIseq_Data_RotateVectors.m
    └── SetUp/
        ├── startup.m
        └── sms-recon/
            ├── initializePaths.m
            └── ...
```

---

## Requirements

This pipeline is intended to run in MATLAB. The required toolbox and external software may depend on the reconstruction and inversion environment used in your local setup.

Typical requirements include:

* MATLAB
* Existing MRE reconstruction utilities in `sms-recon`
* FSL, if registration or NIfTI-based processing steps require FSL commands
* Access to raw spiral or EPI MRE data
* MRE inversion software or downstream pipeline for stiffness reconstruction

Before running the scripts, make sure all required paths are initialized using the provided `startup.m` and `initializePaths.m` files.

---

# Spiral MRE Pipeline

The spiral pipeline performs navigator-based motion estimation, motion-corrected multiband reconstruction, and wave field correction of the reconstructed displacement maps.

## Step 1. Initialize the Spiral Environment

In MATLAB, navigate to the `Spiral` folder and run:

```matlab
startup
```

Then navigate to `Spiral/sms-recon` and run:

```matlab
initializePaths
```

These scripts add the required reconstruction, utility, and processing folders to the MATLAB path.

---

## Step 2. Reconstruct Navigator Images

Run the following script under:

```text
Spiral/sms-recon/Routines/
```

Script:

```matlab
prepMultibandMRE
```

This step reconstructs navigator images from the spiral multiband MRE acquisition. These navigator images are later used for motion estimation.

Expected output:

* Navigator image series
* Intermediate reconstruction files needed for registration

---

## Step 3. Estimate Motion from Navigator Images

Run:

```matlab
calc2Dregistrations_ZS
```

Location:

```text
Spiral/sms-recon/Routines/
```

This script performs 2D registration on the navigator images and estimates motion parameters for each relevant repetition, shot, phase offset, or motion state depending on the acquisition structure.

Expected output:

* Transformation matrices
* Registration outputs
* Motion parameter records, such as rotation and translation estimates

These transformation matrices are used by the motion-corrected reconstruction step.

---

## Step 4. Reconstruct Motion-Corrected Spiral MRE Images

Run:

```matlab
reconMultibandMRE_ZS
```

Location:

```text
Spiral/sms-recon/Routines/
```

This script reconstructs the multiband spiral MRE data while incorporating the transformation matrices estimated in Step 3.

Expected output:

* Motion-corrected complex MRE images
* Reconstructed MRE image series prepared for displacement processing

---

## Step 5. Perform Wave Field Correction and Prepare Motion Maps for MRE Inversion

Run the following scripts in order under:

```text
Spiral/sms-recon/Utilities/MRE/
```

### 5.1 Initial MRE Processing and Vector Rotation

```matlab
proc_mbmre_1_ZS_RotateVectors
```

This script performs the first stage of MRE processing and applies vector rotation to account for motion-related changes in displacement vector orientation.

### 5.2 Manual Masking

```matlab
manual_masking
```

This step is used to generate or refine the brain mask required for downstream MRE processing and inversion.

### 5.3 Final MRE Processing

```matlab
proc_mbmre_2_ZS_RotateVectors
```

This script completes the wave field correction and prepares the final motion maps for MRE inversion.

Expected output:

* Wave-field-corrected motion maps
* Masked displacement maps
* Motion-corrected and vector-rotated MRE outputs ready for inversion

---

# EPI MRE Pipeline

The EPI pipeline performs motion correction and wave field correction section by section using a single main processing script.

## Step 1. Initialize the EPI Environment

In MATLAB, navigate to:

```text
EPI/SetUp/
```

Run:

```matlab
startup
```

Then navigate to:

```text
EPI/SetUp/sms-recon/
```

Run:

```matlab
initializePaths
```

These scripts initialize the required paths for EPI MRE processing.

---

## Step 2. Run the EPI Motion and Wave Field Correction Script

Run the following script section by section:

```matlab
Siemens_EPIseq_Data_RotateVectors
```

This script performs the EPI processing workflow, including:

1. Loading EPI MRE data.
2. Applying motion correction.
3. Rotating displacement vectors according to the estimated motion.
4. Generating wave-field-corrected motion maps.
5. Preparing outputs for MRE inversion.

Expected output:

* Motion-corrected EPI MRE images
* Wave-field-corrected displacement maps
* Final MRE motion maps ready for inversion

---

# Summary of Workflow

## Spiral

```text
Spiral/startup.m
        ↓
Spiral/sms-recon/initializePaths.m
        ↓
prepMultibandMRE.m
        ↓
calc2Dregistrations_ZS.m
        ↓
reconMultibandMRE_ZS.m
        ↓
proc_mbmre_1_ZS_RotateVectors.m
        ↓
manual_masking.m
        ↓
proc_mbmre_2_ZS_RotateVectors.m
        ↓
Motion-corrected and wave-field-corrected MRE maps ready for inversion
```

## EPI

```text
EPI/SetUp/startup.m
        ↓
EPI/SetUp/sms-recon/initializePaths.m
        ↓
Siemens_EPIseq_Data_RotateVectors.m
        ↓
Motion-corrected and wave-field-corrected EPI MRE maps ready for inversion
```

---

# Notes

* The spiral workflow uses navigator images to estimate motion before motion-corrected reconstruction.
* The EPI workflow is currently implemented as a section-by-section script.
* The vector rotation step is important because subject motion changes the orientation of the measured displacement vectors relative to the anatomical frame.
* Manual mask inspection is recommended before MRE inversion.
* Output files should be checked visually before proceeding to stiffness inversion.

---

# Suggested Quality Control

After running the correction pipeline, visually inspect:

* Navigator registration quality
* Estimated transformation matrices
* Motion-corrected magnitude images
* Phase/displacement maps before and after correction
* Brain mask quality
* Final wave-field-corrected motion maps

Recommended tools include MATLAB visualization, FSLeyes, or other NIfTI-compatible image viewers depending on the output file format.

---

# Citation

If you use this pipeline, please cite the associated COMET motion correction work.

```text
Citation information to be added.
```

---

# Contact

For questions, issues, or collaboration inquiries, please contact the repository maintainer.
