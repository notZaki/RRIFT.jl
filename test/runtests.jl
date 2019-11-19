using RRIFT
using Test

@testset "RRIFT" begin
    (vfa_folders, dce_folders) = download_invivo_studies(study=RRIFT.gbm_study_uids[1], destination_folder = "./data/tcga-gbm")
    mask_folder = download_invivo_masks(destination_folder = "./data/tcga-gbm-masks")

    relaxation = compute_concentration(vfa_folder = vfa_folders[1], dce_folder=dce_folders[1])
    masks = get_mask(study=RRIFT.gbm_study_uids[1], mask_folder=mask_folder)

    Ct = apply_mask(data=relaxation.ct, mask = masks.tumour)
    @test size(Ct) == (12895, 70)
end
