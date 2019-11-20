using RRIFT
using Test

@testset "Preprocessing" begin
    (vfa_folders, dce_folders) = download_invivo_studies(RRIFT.gbm_study_uids[1], destination = "./data/tcga-gbm")
    mask_folder = download_invivo_masks(destination = "./data/tcga-gbm-masks")

    relaxation = compute_concentration(vfa_folder = vfa_folders[1], dce_folder=dce_folders[1])
    masks = get_mask(study=RRIFT.gbm_study_uids[1], mask_folder=mask_folder)

    Ct = apply_mask(data=relaxation.ct, mask = masks.tumour)
    @test size(Ct) == (12895, 70)
end

@testset "In-vivo" begin
    mat_folder = "./data/tcga-gbm-mat"
    download_invivo_preprocessed(destination = mat_folder)

    fits = process_patients(mat_folder)

    ccc_kt = ccc(vcat(fits[:kt][:tofts]...), vcat(fits[:kt][:rrift]...), lim=(0, 0.2))
    ccc_ve = ccc(vcat(fits[:ve][:tofts]...), vcat(fits[:ve][:rrift]...), lim=(0, 0.5))
    ccc_vp = ccc(vcat(fits[:vp][:tofts]...), vcat(fits[:vp][:rrift]...), lim=(0, 0.05))
    
    ccc_ktrr = ccc(fits[:kt_rr][:tofts], fits[:kt_rr][:rrift])
    ccc_verr = ccc(fits[:ve_rr][:tofts], fits[:ve_rr][:rrift])
    ccc_keprr = ccc(fits[:kep_rr][:tofts], fits[:kep_rr][:rrift])

    println("CCC [kt, ve vp]:")
    println([ccc_kt, ccc_ve, ccc_vp])
    println("CCC [kt_rr, ve_rr, kep_rr]:")
    println([ccc_ktrr, ccc_verr, ccc_keprr])
end
