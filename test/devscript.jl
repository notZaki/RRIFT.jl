(vfa_folders, dce_folders) = download_invivo_studies(destination = "./data/tcga-gbm")
mask_folder = download_invivo_masks(destination = "./data/tcga-gbm-masks")
mat_files = preprocess_dicom_to_mat(destination = "./data/tcga-gbm-mat", vfa_folders = vfa_folders, dce_folders = dce_folders, mask_folder = mask_folder)

m = process_patients(mat_files)

data = load_preprocessed_mat(mat_files[7])
@. data.ct[data.ct < 0] = 0
@. data.crr[data.crr < 0] = 0
@. data.cp[data.cp < 0] = 0

tmp = fit_model(:extendedtofts, :lls, t=data.t, ct=data.ct, cp=data.cp).estimates
tmprr = fit_cerrm_with_rrift(t=data.t, ct=data.ct, cp=data.cp, crr=data.crr, tail_start=33)

tmprr = fit_cerrm(t=data.t, ct=data.ct, crr=data.crr)
findfirst(data.t.>3)
rrift(t=data.t, cp=data.cp, crr=data.crr, tail_start=35, kep_rr=tmprr.kep_rr)

kt, kep, ve = (fill(NaN, size(data.masks["aif"])) for _=1:3)
kt[data.masks["tumour"]] = tmp.kt
kep[data.masks["tumour"]] = tmp.kep
ve[data.masks["tumour"]] = tmp.ve
kt = autocrop(kt, mask=data.masks["tumour"])
ktr, kepr, ver = (fill(NaN, size(data.masks["aif"])) for _=1:3)
ktr[data.masks["tumour"]] = tmprr.kt
kepr[data.masks["tumour"]] = tmprr.kep
ver[data.masks["tumour"]] = tmprr.ve

heatmap(kt[:,:,4], clim=(0, 0.15), c=:cinferno, yflip=true, aspect_ratio=:equal)

histogram2d(tmp.kt, tmprr.kt, nbins=20, color=cgrad(:cinferno,scale=:log))
