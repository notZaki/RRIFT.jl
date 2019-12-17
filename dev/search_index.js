var documenterSearchIndex = {"docs":
[{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"using RRIFT, Plots, Perfusion, Statistics\n\nENV[\"GKSwstype\"] = \"100\"\ngr()\n\nprintln(\"Fitting\")\nfigopts = (framestyle = :grid, gridalpha=0.5, gridstyle=:dot, linewidth = 2.5, \n        tickfontsize = 11, fg_legend = :transparent, legendfontsize = 11, legend=:bottomright)\n\nlineopts(label::String) = (label = label, xlabel = \"Time [min]\", ylabel = \"[Gd] [mM]\", figopts...)","category":"page"},{"location":"guide/fitting/#In-vivo-analysis-on-glioblastoma-multiforme-(GBM)-1","page":"Fitting","title":"In-vivo analysis on glioblastoma multiforme (GBM)","text":"","category":"section"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"This section will apply the reference region and input function tail (RRIFT) method on the in-vivo data. The Tofts model will also be fitted.","category":"page"},{"location":"guide/fitting/#Downloading-the-pre-processed-data-1","page":"Fitting","title":"Downloading the pre-processed data","text":"","category":"section"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"For simplicity, we will use data that has already been pre-processed and which can be downloaded by:","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"mat_dir = \"./data/tcga-gbm-mat\"\ndownload_invivo_preprocessed(destination = mat_dir)\nmat_files = joinpath.(mat_dir, readdir(mat_dir))","category":"page"},{"location":"guide/fitting/#Single-patient-example-1","page":"Fitting","title":"Single patient example","text":"","category":"section"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"The next few sections will fit the extended Tofts model and the RRIFT method on a single DCE-MRI study.","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"chosen_mat_file = mat_files[7]\nmat_data = load_preprocessed_mat(chosen_mat_file)\n@extract (t, ct, crr, cp, masks) mat_data\n\nprintln(\"\"\"\nLoaded the following variables:\n    - `t` is an $(typeof(t)) with length $(length(t))\n    - `ct` is an $(typeof(ct)) with size $(size(ct))\n    - `crr` is an $(typeof(crr)) with length $(length(crr))\n    - `cp` is an $(typeof(cp)) with length $(length(cp))\n    - `masks` is a $(typeof(masks)) with keys $(keys(masks))\n\"\"\")","category":"page"},{"location":"guide/fitting/#Extended-Tofts-model-fit-1","page":"Fitting","title":"Extended Tofts model fit","text":"","category":"section"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"The extended Tofts model has the form $ Ct(t) = K^{trans} \\cdot Cp(t) \\ast \\exp(-k{ep} \\cdot t) + vp \\cdot Cp(t) $ where Ct$ is the concentration in tissue, C_p is the input function, and t is the time, i.e. ct, cp and, t in the code, respectively. The ast is a convolution while the fitting parameters are K^trans, k_ep, and v_p, along with a derived parameter v_e = K^trans  k_ep (not shown in equation).","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"The input function along with concentration-time data for a single voxel are shown below:","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"plot(t, cp; lineopts(\"Input function\")..., c = :red, legend = :topright)\nsavefig(\"cp.png\"); nothing # hide","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"(Image: cp)","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"The concentration-time data in a single tumour voxel is shown below:","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"single_ct = ct[100, :]\nscatter(t, single_ct; lineopts(\"Concentration-time data in a tumour voxel\")...)\nsavefig(\"single_ct.png\"); nothing # hide","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"(Image: single_ct)","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"Fitting the extended tofts model to the single-voxel curve results in:","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"est_tofts = fit_model(:extendedtofts, :lls, ct = single_ct, t = t, cp = cp).estimates\nfitted_curve = model_tofts(t = t, cp = cp, \n    parameters = (kt = est_tofts.kt[1], kep = est_tofts.kep[1], vp = est_tofts.vp[1]))\n\nscatter(t, single_ct; lineopts(\"Measured curve in single voxel\")...)\nplot!(t, fitted_curve; title = \"Extended Tofts model fit\", lineopts(\"Extended Tofts model fit\")...)\nsavefig(\"tofts_fit.png\"); nothing # hide","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"(Image: tofts_fit) where the fitting parameters are:","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"est_tofts","category":"page"},{"location":"guide/fitting/#Extended-reference-region-model-1","page":"Fitting","title":"Extended reference region model","text":"","category":"section"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"One of the issues with the Tofts model is that it require knowledge of the input function cp. The input function has a sharp initial peak and a fast temporal resolution is needed to accurately measure it.  This requires sacrificing SNR, spatial resolution, and volume coverage, all of which are precious. ","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"An alternative is the reference region model which uses a healthy reference tissue crr as a surrogate for cp. The reference tissue curve is:","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"plot(t, crr; lineopts(\"Reference tissue curve\")..., c=:green)\nsavefig(\"crr.png\"); nothing # hide","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"(Image: crr)","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"Fitting the extended reference region model to the single voxel results in:","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"est_errm = fit_errm(t = t, ct = single_ct, crr = crr)","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"The extended reference region model provides estimates for: K^transK^trans_RR, v_ev_eRR, v_pK^trans_RR,  k_ep, and k_epRR. These are rel_kt, rel_ve, rel_vp, kep and kep_rr in the code, respectively.  The parameters with the RR subscript represent the reference tissue.  In order to get K^trans, v_e, and v_p, we need to know the reference tissue's K^trans_RR and v_eRR.","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"Let's first look at the bright side: the reference region model provides en estimate for k_ep without needing an input function. Unfortunately, this value does not agree with the estimate we got earlier from the Tofts model:","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"println(\"\"\"\nkep estimated with: \n    - Tofts model: $(est_tofts.kep[1])\n    - Ref.Region model: $(est_errm.kep[1])\n\"\"\")","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"Well, let's look at the other bright side: the reference region model provides an estimate for k_epRR.  Each fit on a tumour voxel produces an estimate for k_epRR.  In theory, all fits should estimate the same k_epRR because this parameter described the reference tissue and all fits use the same reference tissue curve. In practice, the estimated value varies due to noise and other fitting artifacts—for example, the voxel's estimated k_epRR is -0.01 which is unphysical. However, most of the fits should be centered around the same value.","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"# Fit all tumour voxels with the extended reference region model\nest_errm_allvoxels = fit_errm(t=t, ct=ct, crr=crr)\n# Plot a histogram of the kep_rr estimates\nhistogram(est_errm_allvoxels.kep_rr[0 .< est_errm_allvoxels.kep_rr .< 2], bins=100, linealpha=0; \n    lineopts(\"Estimated kep_rr\")..., xlabel=\"Estimate kep_rr [1/min]\", ylabel=\"Counts\")\nsavefig(\"hist_keprr.png\"); nothing # hide","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"(Image: hist_keprr)","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"There is a peak in the histogram close to 0.3~0.4. We can estimate a single k_epRR value by considering only the fits with positive estimates and then taking the interquartile mean of k_epRR from those fits.","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"positive_mask = positive_only_mask(est_errm_allvoxels)\nest_kep_rr = interquartile_mean(est_errm_allvoxels.kep_rr[positive_mask])\nprintln(\"Estimated kep_rr: $est_kep_rr\")","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"Now that we have a single estimate for kep_rr, we can re-fit the reference region model but this time we force the fits to have the same kep_rr value. This two-fit approach is called the constrained extended reference region model, and it leads to better agreement with the Tofts fit:","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"# Refit using a fixed kep_rr\nest_cerrm = fit_cerrm(t = t, ct = single_ct, crr = crr, kep_rr = est_kep_rr)\n\nprintln(\"\"\"\nkep estimated with: \n    - Tofts model: $(est_tofts.kep[1])\n    - Constrained Ref.Region model: $(est_cerrm.kep[1])\n\"\"\")","category":"page"},{"location":"guide/fitting/#Reference-region-and-input-function-tail-(RRIFT)-method-1","page":"Fitting","title":"Reference region and input function tail (RRIFT) method","text":"","category":"section"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"Let's return our attention to the fact that the reference region model gives us relative parameters.  In order to get absolute parameters, we need to know the reference tissue parameters: K^trans_RR and v_eRR. This is typically done by using literature-based values of muscle, but K^trans_RR varies substantially between patients and between muscles.","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"The paper proposes RRIFT which takes advantages of two features:","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"We know k_epRR already. This is useful because k_epRR = K^trans_RRv_eRR, so we only need to know either K^trans_RR or v_eRR.\nThe peak part of the input function cp is hard to measure, but the rest of the input function is fairly \"flat\" and could be measured with a slow scan. This \"input function tail\" is plotted next:","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"tail_start = findfirst(t .> 2)\nplot(t, cp; lineopts(\"Input function\")...)\nplot!(t[tail_start:end], cp[tail_start:end]; lineopts(\"Input function tail\")...)\nsavefig(\"tail.png\"); nothing # hide","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"(Image: tail)","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"The equation to estimate K^trans_RR is:","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"K^trans_RR = fracC_RR(t) - C_RR(t_start) + k_epRR cdot int_t_start^t C_RR(tau) dtauint_t_start^t C_p(tau) dtau","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"where t_start is the start of the AIF tail (2 minutes in the above example), t is any timepoint after t_start, and we already estimated k_epRR earlier from the reference region model.","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"The equation can be solved by linear regressing with t = t_start+1  t_end, as shown next:","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"using NumericalIntegration: cumul_integrate\n\ntail_start = findfirst(t .> 2)\ncrr_tail = crr[tail_start:end]\ncp_tail = cp[tail_start:end]\nt_tail = t[tail_start:end]\n\nnumerator = crr_tail .- crr_tail[1] .+ est_kep_rr .* cumul_integrate(t_tail, crr_tail)\ndenominator = cumul_integrate(t_tail, cp_tail)\n\nest_kt_rr = denominator \\ numerator\nest_ve_rr = est_kt_rr / est_kep_rr\nprintln(\"Estimated Ktrans_rr from RRIFT fit: $(round(est_kt_rr, digits=4))\")\nprintln(\"Estimated ve_rr: $(round(est_ve_rr, digits=4))\")\n\nscatter(denominator, numerator; lineopts(\"Data\")..., legend=:bottomright)\nplot!(denominator, denominator .* est_kt_rr; \n    lineopts(\"RRIFT fit\")..., xlabel = \"Denominator\", ylabel=\"Numerator\")\nsavefig(\"rrift_fit.png\"); nothing # hide","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"(Image: rrift_fit)","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"Now we can use the estimated K^trans_RR and v_eRR to convert the relative estimates from the reference region model into absolute estimates.","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"est_tofts = fit_model(:extendedtofts, :lls, t = t, ct = single_ct, cp = cp).estimates\nest_cerrm = fit_cerrm(t = t, ct = single_ct, crr = crr, kep_rr = est_kep_rr)\n\nprintln(\"\"\"\nComparison between Toft and RRM+RRIFT\n    - Tofts Ktrans: $(est_tofts.kt)\n    - RRIFT Ktrans: $(est_cerrm.rel_kt .* est_kt_rr)\n    --------------------------------------------\n    - Tofts ve: $(est_tofts.ve)\n    - RRIFT ve: $(est_cerrm.rel_ve .* est_ve_rr)\n    -------------------------------------------\n    - Tofts vp: $(est_tofts.vp)\n    - RRIFT vp: $(est_cerrm.rel_vp .* est_kt_rr)\nThey're quite similar! :)\n\"\"\")","category":"page"},{"location":"guide/fitting/#Voxel-wise-fitting-1","page":"Fitting","title":"Voxel-wise fitting","text":"","category":"section"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"The above was an example for a single tumour voxel.  This section will use voxel-wise fitting to show that the maps using RRIFT and similar to the Tofts model.","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"tail_start = findfirst(t .> 2)\nest_tofts = fit_model(:extendedtofts, :lls, t = t, ct = ct, cp = cp).estimates\nest_rrift = fit_cerrm_with_rrift(t = t, ct = ct, cp = cp, crr = crr, tail_start = tail_start)\n# Note: I didn't pass `kep_rr` as an input argument above. That's because the function will compute it on its own.\n\nestimates = (tofts = est_tofts, rrift = est_rrift)\n\nmaps = Dict()\nfor param in keys(est_tofts)\n    inner_dict = Dict()\n    for method in (:tofts, :rrift)\n        i_am_the_map = zeros(size(masks[\"tumour\"]))\n        i_am_the_map[masks[\"tumour\"]] .= estimates[method][param]\n        i_am_the_map = RRIFT.crop(i_am_the_map)\n        inner_dict[method] = i_am_the_map\n    end\n    maps[param] = inner_dict\nend\n\nslice = 6\np1 = heatmap(maps[:kt][:tofts][:,:,slice], c=:cinferno, yflip=true, aspect_ratio=:equal, clim=(0, 0.2); lineopts(\"\")..., title=\"Tofts\", axis=nothing, xlabel=\"\", ylabel=\"\")\np2 = heatmap(maps[:kt][:rrift][:,:,slice], c=:cinferno, yflip=true, aspect_ratio=:equal, clim=(0, 0.2); lineopts(\"\")..., title=\"RRIFT\", axis=nothing, xlabel=\"\", ylabel=\"\")\nplot(p1, p2, layout=(1,2))\nsavefig(\"voxelwise.png\"); nothing # hide","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"(Image: voxelwise)","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"println(\"\"\"\nConcordance correlation coefficients between Tofts and RRIFT fits\n    - for ktrans: $(ccc(est_tofts[:kt], est_rrift[:kt], lim=(0, 0.2)))\n    - for ve: $(ccc(est_tofts[:ve], est_rrift[:ve], lim=(0, 0.5)))\n    - for vp: $(ccc(est_tofts[:vp], est_rrift[:vp], lim=(0, 0.05)))\n\"\"\")","category":"page"},{"location":"guide/fitting/#","page":"Fitting","title":"Fitting","text":"println(\"End Fitting\")","category":"page"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"using RRIFT\n\nprintln(\"Downloading\")\nchosen_study_uid = RRIFT.gbm_study_uids[8]\ndicom_folders = download_invivo_studies(chosen_study_uid, destination = \"./data/tcga-gbm-dicom\")","category":"page"},{"location":"guide/downloading/#Downloading-DICOM-files-1","page":"Downloading","title":"Downloading DICOM files","text":"","category":"section"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"The in-vivo evaluation uses publicly available data in The Cancer Genome Atlas - Glioblastoma Multiforme (TCGA-GBM) collection from The Cancer Imaging Archive (TCIA).","category":"page"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"The TCGA-GBM collection contains nearly 600 studies with over 5,000 imaging series, however not all of these contain DCE-MRI scans.  The RRIFT manuscript used 8 DCE-MRI studies and their unique identifiers (Study Instance UIDs) were included in supplementary materials table S1. These UIDs can be used to download the DICOM files.","category":"page"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"RRIFT.gbm_study_uids","category":"page"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"Each study contains multiple imaging series. For our purposes, we need the DCE-MRI series along with the variable flip angle (VFA) series.  The respective DICOM files for all 8 studies can be downloaded into a destination folder by","category":"page"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"dicom_folders = download_invivo_studies(destination = \"./data/tcga-gbm\")","category":"page"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"The above function will automatically identify the DCE-MRI and VFA series for each study.","category":"page"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"To only download a single study, pass a study UID to the function.","category":"page"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"chosen_study_uid = RRIFT.gbm_study_uids[8]\ndicom_folders = download_invivo_studies(chosen_study_uid, \n    destination = \"./data/tcga-gbm-dicom\")","category":"page"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"note: Note\nThe function will not download anything if the destination folder already contains the dicom files.  To force a download, pass overwrite = true as an argument to the function. This tip applies to all other download functions as well.","category":"page"},{"location":"guide/downloading/#Downloading-masks-1","page":"Downloading","title":"Downloading masks","text":"","category":"section"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"Masks/contours for the tissue of interest (tumour), temporalis muscle (reference tissue), and arterial input function were manually drawn and saved as .mat files. They can be downloaded by:","category":"page"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"download_invivo_masks(destination = \"./data/tcga-gbm-masks\")","category":"page"},{"location":"guide/downloading/#Downloading-Pre-processed-.mat-data-1","page":"Downloading","title":"Downloading Pre-processed .mat data","text":"","category":"section"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"The input data required a couple of pre-processing steps which include:","category":"page"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"Loading the VFA and DCE-MRI DICOM files\nComputing T1 maps from the VFA signal\nConverting DCE-MRI signal to tracer concentration\nApplying the masks to extract curves for the tumour, reference tissue, and input function.","category":"page"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"To save time, the product of these post-processing steps has been saved in .mat files which can be downloaded by:","category":"page"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"download_invivo_preprocessed(destination = \"./data/tcga-gbm-mat\")","category":"page"},{"location":"guide/downloading/#","page":"Downloading","title":"Downloading","text":"println(\"End Downloading\")","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"using RRIFT, Perfusion, Statistics, Plots\n\nENV[\"GKSwstype\"] = \"100\"\ngr()\n\nprintln(\"Preprocessing\")","category":"page"},{"location":"guide/preprocessing/#Pre-processing-1","page":"Pre-processing","title":"Pre-processing","text":"","category":"section"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"Before applying RRIFT, there's a couple of pre-processing steps including:","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"Load imaging data from the VFA DICOM files & compute T1 maps\nLoad imaging data form the DCE DICOM files & convert the DCE-MRI signal into tracer concentration\nExtract signal-time curves from the tumour, muscle, and artery for subsequent model fitting","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"Download DICOM files","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"The pre-processing steps will be shown for a single patient. First, the DICOM files must be downloaded.","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"chosen_study_uid = RRIFT.gbm_study_uids[8]\n\ndicom_folders = download_invivo_studies(chosen_study_uid, \n    destination = \"./data/tcga-gbm-dicom\")\n\n# Extract the vfa and dce folders from `dicom_folders`\nvfa_folder = dicom_folders.vfa_folders[1]\ndce_folder = dicom_folders.dce_folders[1]\n\nprintln(\"VFA: $vfa_folder\")\nprintln(\"DCE: $dce_folder\")","category":"page"},{"location":"guide/preprocessing/#T1-mapping-1","page":"Pre-processing","title":"T1 mapping","text":"","category":"section"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"T1 mapping is provided by the fit_relaxation function.  There are three algorithms for T1 mapping:","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"Non-linear least squares fitting of the spoiled gradient echo equation\nLinear least squares fitting with DESPOT1\nIterative fitting with NOVIFAST","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"The fitting algorithm is selected by passing either :nls, :despot, or :novifast as the first argument. This example uses DESPOT1 because that is what the paper used. I wasn't aware of NOVIFAST when I wrote the paper, or else I probably would've used it instead of DESPOT1.","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"# Load VFA data\nvfa = RRIFT.load_vfa_dicom(folder = vfa_folder)\n\nprintln(\"\"\"\nSome information for VFA data:\n    - It is a named tuple with keys: $(keys(vfa))\n    - Number of flip angles: $(length(vfa.angles))\n    - Value of flip angles, in degrees: $(sort(round.(rad2deg.(vfa.angles))))\n    - Size of signal data: $(size(vfa.signal))\n    - Repetition time, in ms: $(vfa.TR)\n\"\"\")\n\n# Compute T1 maps using DESPOT1\nrelaxation_maps = fit_relaxation(:despot; vfa...).estimates\n\n# The variable `relaxation_maps` contains the keys: T1 & M0\nkeys(relaxation_maps)","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"note: Note\nThe fit_relaxation function expects the first argument to be a symbol followed by keyword arguments.  The keywords are signal, angles, and TR.  In others words, it could've been written as: fit_relaxation(:despot, signal = vfa.signal, angles = vfa.angles, TR = vfa.TR). However, since the keys in vfa match the function's keywords, we can just splat it with .... ","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"# Remove empty rows/columns for display purposes\nM0 = crop(relaxation_maps.M0)\nT1 = crop(relaxation_maps.T1)\n\nslice = 8\nfigopts = (c = :cinferno, yflip = true, aspect_ratio = :equal, axis = nothing)\np1 = heatmap(M0[:,:,slice]; clim=(0, 15000), title=\"M0\", figopts...)\np2 = heatmap(T1[:,:,slice]; clim=(500, 3000), title=\"T1\", figopts...)\nplot(p1, p2, layout=(1,2))\nsavefig(\"relaxation_maps.png\"); nothing # hide","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"(Image: relaxation_maps)","category":"page"},{"location":"guide/preprocessing/#DCE-MRI-signal-to-concentration-conversion-1","page":"Pre-processing","title":"DCE-MRI signal to concentration conversion","text":"","category":"section"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"The T1 map is used to convert the DCE-MRI signal into tracer concentration.","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"# Load the DCE-MRI data\ndce = RRIFT.load_dce_dicom(folder = dce_folder)\n\nprintln(\"\"\"\nSome information for DCE data:\n    - It is a named tuple with keys: $(keys(dce))\n    - Size of signal data: $(size(dce.signal))\n    - Number of timepoints/frames: $(length(dce.timepoints))\n    - Flip angle, in degrees: $(rad2deg(dce.angle))\n    - Repetition time, in ms: $(dce.TR)\n\"\"\")\n\nr1 = 3.3/1000 # Relaxivity of Gd-DTPA at 3T, in mM/ms. Ref: PMID 16481903\nBAF = 3 # Bolus arrival frame, i.e. the frame at which the tracer arrives in the imaging volume\nR10 = 1 ./ relaxation_maps.T1 # The function wants R1 instead of T1\n\nconcentration = signal_to_concentration(dce.signal; \n    angle = dce.angle, TR = dce.TR, R10 = R10, BAF = BAF, r1 = r1)\nprintln(\"Concentration size: $(size(concentration))\")","category":"page"},{"location":"guide/preprocessing/#Wrapper-function-for-T1-mapping-and-signal-concentration-conversion-1","page":"Pre-processing","title":"Wrapper function for T1 mapping and signal-concentration conversion","text":"","category":"section"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"The functions in the previous section were shown for educational purposes.  The compute_concentration function wraps the previous steps together:","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"computed = compute_concentration(vfa_folder = vfa_folder, dce_folder = dce_folder)\n\nprintln(\"\"\"\nSome information for the computed maps:\n    - It is a named tuple with keys: $keys(computed)\n    - Does it have the same T1 map as we computed earlier? $(all(@. (computed.T1 == relaxation_maps.T1)[!isnan(computed.T1)]))\n    - Does it have the same concentration as computed earlier? $(all(@. (computed.ct == concentration)[!isnan(concentration)]))\n\"\"\")","category":"page"},{"location":"guide/preprocessing/#Applying-masks-1","page":"Pre-processing","title":"Applying masks","text":"","category":"section"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"Pharmacoknetic modelling of DCE-MRI requires concentration-time curves from:","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"the tissue of interest (tumour in our cases) denoted as ct\nthe feeding artery (arterial input function / AIF) denoted as cp\ntechnically, it should be cb and we'll actually be using a vein instead of artery\na healthy reference tissue (muscle) denoted as crr","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"Not all models require the same curves. For example, the well-established Tofts model only needs ct and cp, whereas the reference region model needs ct and crr.","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"Contours/masks for the tumour/muscle/blood-vessel were manually drawn and they can be downloaded by:","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"mask_folder = \"./data/tcga-gbm-masks\"\ndownload_invivo_masks(destination = mask_folder)","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"Masks can be loaded and applied by:","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"# Get masks for the chosen study\nmask = get_mask(study = chosen_study_uid, mask_folder = mask_folder)\nprintln(\"Mask has keys: $(keys(mask))\")\n\nhematocrit = 0.4 # For AIF (assumed value)\n\nt = computed.t\nct = apply_mask(data = computed.ct, mask = mask.tumour)\ncp = apply_mask(data = computed.ct, mask = mask.aif) ./ (1 - 0.4)\ncrr = apply_mask(data = computed.ct, mask = mask.muscle)\n\n# cp and crr are averaged because only a single representative curve is needed from each\ncp = vec(mean(cp, dims=1))\ncrr = vec(mean(crr, dims=1))\nnothing # Suppressing output","category":"page"},{"location":"guide/preprocessing/#(Another)-Wrapper-function-for-all-pre-processing-steps-1","page":"Pre-processing","title":"(Another) Wrapper function for all pre-processing steps","text":"","category":"section"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"The pre-processing steps—i.e. T1 mapping, signal-concentration conversion, masking—are needed for each study. Rather than repeating these steps every time, a wrapper function named preprocess_dicom_to_mat is used.  This function applies all of the preceding steps and saves the result as a MATLAB-compatible .mat file.","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"preprocessed_mat_files = preprocess_dicom_to_mat(\n    destination = \"./data/tcga-gbm-mat-test\", \n    dicom_folders = dicom_folders, \n    mask_folder = mask_folder)","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"The saved .mat file can be loaded back into Julia by","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"mat_data = load_preprocessed_mat(preprocessed_mat_files[1])\n\nprintln(\"\"\"\nSome information for the loaded mat_data\n    - It is a named tuple with keys: $(keys(mat_data))\n    - `masks` is a dictionary with keys: $(keys(mat_data.masks))\n    - `relaxation` is a dictionary with keys: $(keys(mat_data.relaxation))\n\"\"\")","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"All 8 DCE-MRI studies can be pre-processed by running the following three lines:","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"dicom_folders = download_invivo_studies(destination = \"./data/tcga-gbm-dicom\")\nmask_folder = download_invivo_masks(destination = \"./data/tcga-gbm-masks\")\npreprocess_dicom_to_mat(destination = \"./data/tcga-gbm-mat\", \n    dicom_folders = dicom_folders, mask_folder = mask_folder)","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"To save time, these steps were already done and the pre-processed files can be downloaded directly by:","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"mat_dir = \"./data/tcga-gbm-mat\"\ndownload_invivo_preprocessed(destination = mat_dir)\n\n# Confirm that files were downloaded:\nmat_files = readdir(mat_dir)","category":"page"},{"location":"guide/preprocessing/#","page":"Pre-processing","title":"Pre-processing","text":"println(\"End Preprocessing\")","category":"page"},{"location":"#Introduction-1","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"This is a Julia module which implements the reference region and input function tail (RRIFT) method described in:","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"Zaki Ahmed & Ives R. Levesque. Pharmacokinetic modeling of dynamic contrast‐enhanced MRI using a reference region and input function tail. Magn Reson Med. 2020; 83: 286– 298. https://doi.org/10.1002/mrm.27913","category":"page"},{"location":"#Installation-1","page":"Introduction","title":"Installation","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"The module can be installed by","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"julia> ]add https://github.com/notZaki/RRIFT.jl.git#master","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"Once installed it can be used by","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"julia> using RRIFT","category":"page"},{"location":"#Usage-1","page":"Introduction","title":"Usage","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"The next sections will go over in-vivo evaluation of RRIFT on in-vivo data. An interactive version of these sections is available on mybinder by clicking on the badge below:  ","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"(Image: Binder)","category":"page"}]
}
