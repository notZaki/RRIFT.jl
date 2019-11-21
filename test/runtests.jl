using RRIFT
using Test

@testset "in-vivo" begin
    test_study = RRIFT.gbm_study_uids[1]
    # Run download functions twice to check for overwrite protection
    dicom_folders = download_invivo_studies(test_study, destination = "./data/tcga-gbm")
    download_invivo_studies(test_study, destination = "./data/tcga-gbm")
    mask_folder = download_invivo_masks(destination = "./data/tcga-gbm-masks")
    download_invivo_masks(destination = "./data/tcga-gbm-masks")

    mat_file_test = preprocess_dicom_to_mat(destination = "./data/tcga-gbm-test", dicom_folders = dicom_folders, mask_folder = mask_folder)

    mat_folder = "./data/tcga-gbm-mat"
    download_invivo_preprocessed(destination = mat_folder)
    
    test_fits_a = process_patients(mat_file_test)
    test_fits_b = process_patients(joinpath(mat_folder, test_study * ".mat"))

    thresh = 1e-9
    for param in (:kt, :kep, :ve, :vp, :kt_rr, :kep_rr, :ve_rr), method in (:tofts, :rrift)
        Δ = maximum(abs.(test_fits_a[param][method][1] .- test_fits_b[param][method][1]))
        println(string(param) * " " * string(method) * ": " * string(Δ))
        @test thresh > Δ 
    end

    fits = process_patients(mat_folder)

    ccc_kt = ccc(fits[:kt], lim=(0, 0.2))
    ccc_ve = ccc(fits[:ve], lim=(0, 0.5))
    ccc_vp = ccc(fits[:vp], lim=(0, 0.05))
    
    ccc_ktrr = ccc(fits[:kt_rr])
    ccc_verr = ccc(fits[:ve_rr])
    ccc_keprr = ccc(fits[:kep_rr])

    println("CCC [kt, ve vp]:")
    println([ccc_kt[:rrift][:tofts], ccc_ve[:rrift][:tofts], ccc_vp[:rrift][:tofts]])
    println("CCC [kt_rr, ve_rr, kep_rr]:")
    println([ccc_ktrr[:rrift][:tofts], ccc_verr[:rrift][:tofts], ccc_keprr[:rrift][:tofts]])
end
