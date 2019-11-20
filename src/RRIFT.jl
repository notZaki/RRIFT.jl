module RRIFT

using Statistics

include("utils.jl")
export apply_mask, autocrop, negatives_to_zero!, ccc

using CancerImagingArchive: series, images
using DICOM: dcm_parse, lookup
include("download_gbm_data.jl")
export download_invivo_study, download_invivo_studies, download_invivo_masks, download_invivo_preprocessed
export get_mask

using MAT: matread, matwrite
using Perfusion: fit_relaxation, signal_to_concentration 
include("preprocess.jl")
export compute_concentration, preprocess_dicom_to_mat, load_preprocessed_mat

using NumericalIntegration: cumul_integrate, TrapezoidalFast
using Perfusion: interquartile_mean, resolve_mask_size, positive_only_mask
include("cerrm.jl")
export fit_cerrm
using Statistics: mean, quantile, std

using Perfusion: fit_model
include("analysis.jl")
export process_patients

function fit_rrift(; crr::AbstractVector, cp::AbstractVector, t::AbstractVector, tail_start::Int, kep_rr::Number)
    @assert length(crr) == length(cp) == length(t)
    @assert tail_start < length(crr)
    crr_tail = crr[tail_start:end]
    cp_tail = cp[tail_start:end]
    t_tail = t[tail_start:end]
    numerator = crr_tail .- crr_tail[1] .+ kep_rr .* cumul_integrate(t_tail, crr_tail, TrapezoidalFast())
    denominator = cumul_integrate(t_tail, cp_tail, TrapezoidalFast())
    kt_rr = denominator \ numerator
    return kt_rr
end
function fit_cerrm_with_rrift(; crr, cp, t, ct, tail_start)
    cerrm = fit_cerrm(crr = crr, ct = ct, t = t)
    kep_rr = cerrm.kep_rr
    kt_rr = fit_rrift(t = t, cp = cp, crr = crr, kep_rr = kep_rr, tail_start = tail_start)
    ve_rr = kt_rr / kep_rr
    
    kt = @. cerrm.rel_kt * kt_rr
    ve = @. cerrm.rel_ve * ve_rr
    vp = @. cerrm.rel_vp * kt_rr
    kep = cerrm.kep

    return (kt = kt, kep = kep, ve = ve, vp = vp, kep_rr = kep_rr, kt_rr = kt_rr, ve_rr = ve_rr)
end
export fit_rrift, fit_cerrm_with_rrift
end # module
