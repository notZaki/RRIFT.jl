module RRIFT

using Statistics

include("utils.jl")
export apply_mask, crop, negatives_to_zero!, ccc

using CancerImagingArchive: series, images
using DICOM: dcm_parse, lookup
using Perfusion: make_folder
include("download_gbm_data.jl")
export download_invivo_studies, download_invivo_masks, download_invivo_preprocessed
export get_mask

using MAT: matread, matwrite
using Perfusion: @extract, fit_relaxation, signal_to_concentration 
include("preprocess.jl")
export compute_concentration, preprocess_dicom_to_mat, load_preprocessed_mat

using Statistics: mean, quantile, std
using Perfusion: fit_model
include("analysis.jl")
export process_patients

import Perfusion

fit_errm(; kwargs...) = Perfusion.fit_errm_lls(; kwargs...).estimates
fit_cerrm(; kwargs...) = Perfusion.fit_cerrm_lls(; kwargs...).estimates

function fit_cerrm_with_rrift(; kwargs...)
    return Perfusion.fit_rrift_with_cerrm(; kwargs...).estimates
end

positive_only_mask = Perfusion.positive_only_mask

export fit_errm, fit_cerrm, fit_cerrm_with_rrift
export positive_only_mask

end # module
