(vfa_folders, dce_folders) = download_invivo_studies(destination = "./data/tcga-gbm")
mask_folder = download_invivo_masks(destination = "./data/tcga-gbm-masks")
mat_files = preprocess_dicom_to_mat(overwrite = true, destination = "./data/tcga-gbm-mat", vfa_folders = vfa_folders, dce_folders = dce_folders, mask_folder = mask_folder)

data = load_preprocessed_mat(mat_files[7])
tmp = fit_model(:tofts, :lls, t=data.t, ct=data.ct, cp=data.cp./(1-0.4)).estimates

kt, kep, ve = (zeros(size(data.masks["aif"])) for _=1:3)
kt[data.masks["tumour"]] = tmp.kt
kep[data.masks["tumour"]] = tmp.kep
ve[data.masks["tumour"]] = tmp.ve
