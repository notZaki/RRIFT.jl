function process_patients(mat_files::AbstractVector{String})
    (kt, kep, ve, vp, kt_rr, kep_rr, ve_rr, vp_rr) = (Dict(:tofts => [], :rrift => []) for _=1:8)
    mask = []
    for file in mat_files
        data = load_preprocessed_mat(file)
        push!(mask, data.masks["tumour"])
        ct = data.ct
        cp = data.cp
        crr = data.crr
        t = data.t
        negatives_to_zero!(ct)
        negatives_to_zero!(cp)
        negatives_to_zero!(crr)

        tail_start = findfirst(t .> 3)

        tofts = fit_model(:extendedtofts, :lls, t = t, ct = ct, cp = cp).estimates
        rrift = fit_cerrm_with_rrift(t = t, ct = ct, cp = cp, crr = crr, tail_start = tail_start)

        push!(kt[:tofts], tofts.kt)
        push!(ve[:tofts], tofts.ve)
        push!(vp[:tofts], tofts.vp)
        push!(kep[:tofts], tofts.kep)
        push!(kt[:rrift], rrift.kt)
        push!(ve[:rrift], rrift.ve)
        push!(vp[:rrift], rrift.vp)
        push!(kep[:rrift], rrift.kep)

        tofts = fit_model(:extendedtofts, :lls, t = t, ct = crr, cp = cp).estimates
        
        push!(kt_rr[:tofts], tofts.kt[1])
        push!(ve_rr[:tofts], tofts.ve[1])
        push!(vp_rr[:tofts], tofts.vp[1])
        push!(kep_rr[:tofts], tofts.kep[1])
        push!(kt_rr[:rrift], rrift.kt_rr[1])
        push!(ve_rr[:rrift], rrift.ve_rr[1])
        push!(vp_rr[:rrift], 0.0)
        push!(kep_rr[:rrift], rrift.kep_rr[1])
    end
    return (kt=kt, kep=kep, ve=ve, vp=vp, kt_rr=kt_rr, kep_rr=kep_rr, ve_rr=ve_rr, vp_rr=vp_rr, mask=mask) 
end

function process_patients(file_or_folder_name::AbstractString)
    if isfile(file_or_folder_name)
        return process_patients([file_or_folder_name])
    elseif isdir(file_or_folder_name)
        return process_patients(joinpath.(file_or_folder_name, readdir(file_or_folder_name)))
    else
        error("$file_or_folder_name is neither a valid file nor directory")
    end
    return
end
