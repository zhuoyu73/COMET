# COMET Motion Correction and Wave Field Reorientation

This repository provides MATLAB pipelines for **COMET motion correction** and **wave field reorientation** for brain Magnetic Resonance Elastography (MRE) data. The workflow supports both **multishot, multiband, jointly-reconstructed MRE** and **single-shot MRE** acquisitions.

The goal of the pipeline is to:

1. Estimate motion from navigator or image-domain data.
2. Incorporate the estimated motion into image reconstruction or post-processing.
3. Rotate displacement vector fields into a consistent anatomical frame.
4. Generate motion-corrected and wave-field-reoriented MRE displacement maps ready for MRE inversion.

---

## Repository Structure

```text
COMET/
├── Iterative Motion Correction/
│   ├── startup.m
│   └── sms-recon/
│       ├── initializePaths.m
│       ├── Routines/
│       │   ├── prep_recon.m
│       │   ├── motion_extraction.m
│       │   └── recon_motion_correction.m
│       └── Utilities/
│           └── MRE/
│               ├── proc_mbmre_1_reorientation.m
│               ├── manual_masking.m
│               └── proc_mbmre_2_reorientation.m
│
└── Post-Reconstruction Correction/
    └── post_recon_motion_correction_reorientation.m
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
* Access to raw multishot or single-shot MRE data
* MRE inversion software or downstream pipeline for stiffness reconstruction

Before running the scripts, make sure all required paths are initialized using the provided `startup.m` and `initializePaths.m` files.

---

# Iterative Reconstruction Correction Pipeline

The Iterative Reconstruction Correction pipeline performs navigator-based motion estimation, motion-corrected multiband reconstruction, and wave field reorientation.

## Step 1. Initialize the Environment

In MATLAB, navigate to the `Iterative Reconstruction Correction` folder and run:

```matlab
startup
```

Then navigate to `Iterative Reconstruction Correction/sms-recon` and run:

```matlab
initializePaths
```

These scripts add the required reconstruction, utility, and processing folders to the MATLAB path.

---

## Step 2. Reconstruct Navigator Images

Run the following script under:

```text
Iterative Reconstruction Correction/sms-recon/Routines/
```

Script:

```matlab
prep_recon
```

This step reconstructs navigator images from the multiband MRE acquisition. These navigator images are later used for motion estimation.

Expected output:

* Navigator image series
* Intermediate reconstruction files needed for registration

---

## Step 3. Estimate Motion from Navigator Images

Run:

```matlab
motion_extraction
```

Location:

```text
Iterative Reconstruction Correction/sms-recon/Routines/
```

This script performs 2D registration on the navigator images and estimates motion parameters for each relevant repetition, shot, phase offset, or motion state depending on the acquisition structure.

Expected output:

* Transformation matrices
* Registration outputs
* Motion parameter records, such as rotation and translation estimates

These transformation matrices are used by the motion-corrected reconstruction step.

---

## Step 4. Reconstruct Motion-Corrected Multi-Shot MRE Images

Run:

```matlab
recon_motion_correction
```

Location:

```text
Iterative Reconstruction Correction/sms-recon/Routines/
```

This script reconstructs the multiband multishot MRE data while incorporating the transformation matrices estimated in Step 3.

Expected output:

* Motion-corrected complex MRE images

---

## Step 5. Perform Wave Field Reorientation and Prepare Motion Maps for MRE Inversion

Run the following scripts in order under:

```text
Iterative Reconstruction Correction/sms-recon/Utilities/MRE/
```

### 5.1 Initial MRE Processing and Vector Rotation

```matlab
proc_mbmre_1_reorientation
```

This script performs the first stage of MRE processing and applies vector rotation to account for motion-related changes in displacement vector orientation.

### 5.2 Manual Masking

```matlab
manual_masking
```

This step is used to generate or refine the brain mask required for downstream MRE processing and inversion.

### 5.3 Final MRE Processing

```matlab
proc_mbmre_2_reorientation
```

This script completes the wave field reorientation and prepares the final motion maps for MRE inversion.

Expected output:

* Manual mask
* Motion-corrected and vector-rotated MRE outputs ready for inversion

---

# Post-Recon Correction Pipeline

The post-recon pipeline performs motion correction and wave field reorientation section by section using a single main processing script.

## Step 1. Initialize the Environment

In MATLAB, navigate to:

```text
Post-Reconstruction Correction/SetUp/
```

Run:

```matlab
startup
```

Then navigate to:

```text
Post-Reconstruction Correction/SetUp/sms-recon/
```

Run:

```matlab
initializePaths
```

These scripts initialize the required paths for Post-Reconstruction Correction.

---

## Step 2. Run the Post-Reconstruction Motion and Wave Field Reorientation Script

Run the following script section by section:

```matlab
single_shot_motion_correction_reorientation
```

This script performs the Post-Reconstruction processing workflow, including:

1. Loading magnitude and phase MRE data & changing the parameters.
2. Applying motion correction.
3. Rotating displacement vectors according to the estimated motion.
4. Generating wave-field-reoriented motion maps.
5. Preparing outputs for MRE inversion.

Expected output:

* Motion-corrected MRE images
* Final wave-field reoriented MRE motion maps ready for inversion

---

# Summary of Workflow

## Iterative Reconstruction Correction

```text
Iterative Reconstruction Correction/startup.m
        ↓
Iterative Reconstruction Correction/sms-recon/initializePaths.m
        ↓
prep_recon.m
        ↓
motion_extraction.m
        ↓
recon_motion_correction.m
        ↓
proc_mbmre_1_reorientation.m
        ↓
manual_masking.m
        ↓
proc_mbmre_2_reorientation.m
        ↓
Motion-corrected and wave-field-reoriented MRE maps ready for inversion
```

## Post-Reconstruction Correction

```text
Post-Reconstruction Correction/SetUp/startup.m
        ↓
Post-Reconstruction Correction/SetUp/sms-recon/initializePaths.m
        ↓
single_shot_motion_correction_reorientation.m
        ↓
Motion-corrected and wave-field-reoriented MRE maps ready for inversion
```

---

# Citation

If you use this pipeline, please cite the associated COMET motion correction work.

```text
Manuscript is Under Review.
```

---

# Contact

For questions, issues, or collaboration inquiries, please contact Zhuoyu Shi at zs2708@columbia.edu.
