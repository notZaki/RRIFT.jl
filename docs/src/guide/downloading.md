```@setup ex
using RRIFT

println("Downloading")
chosen_study_uid = RRIFT.gbm_study_uids[8]
dicom_folders = download_invivo_studies(chosen_study_uid, destination = "./data/tcga-gbm-dicom")
```

# Downloading DICOM files

The in-vivo evaluation uses publicly available data in [The Cancer Genome Atlas - Glioblastoma Multiforme (TCGA-GBM) collection](https://wiki.cancerimagingarchive.net/display/Public/TCGA-GBM) from The Cancer Imaging Archive (TCIA).

The TCGA-GBM collection contains nearly 600 studies with over 5,000 imaging series, however not all of these contain DCE-MRI scans. 
The RRIFT manuscript used 8 DCE-MRI studies and their unique identifiers (Study Instance UIDs) were included in supplementary materials table S1.
These UIDs can be used to download the DICOM files.

```@repl ex
RRIFT.gbm_study_uids
```

Each study contains multiple imaging series.
For our purposes, we need the DCE-MRI series along with the variable flip angle (VFA) series. 
The respective DICOM files for all 8 studies can be downloaded into a `destination` folder by
```julia
dicom_folders = download_invivo_studies(destination = "./data/tcga-gbm")
```
The above function will automatically identify the DCE-MRI and VFA series for each study.

To only download a single study, pass a study UID to the function.
```@repl ex
chosen_study_uid = RRIFT.gbm_study_uids[8]
dicom_folders = download_invivo_studies(chosen_study_uid, 
    destination = "./data/tcga-gbm-dicom")
```

!!! note
    
    The function will not download anything if the `destination` folder already contains the dicom files. 
    To force a download, pass `overwrite = true` as an argument to the function.
    This tip applies to all other download functions as well.

# Downloading masks

Masks/contours for the tissue of interest (tumour), temporalis muscle (reference tissue), and arterial input function were manually drawn and saved as .mat files.
They can be downloaded by:
```@repl ex
download_invivo_masks(destination = "./data/tcga-gbm-masks")
```

# Downloading Pre-processed .mat data

The input data required a couple of pre-processing steps which include:

1. Loading the VFA and DCE-MRI DICOM files
2. Computing T1 maps from the VFA signal
3. Converting DCE-MRI signal to tracer concentration
4. Applying the masks to extract curves for the tumour, reference tissue, and input function.

To save time, the product of these post-processing steps has been saved in .mat files which can be downloaded by:
```@repl ex
download_invivo_preprocessed(destination = "./data/tcga-gbm-mat")
```

```@setup ex
println("End Downloading")
```