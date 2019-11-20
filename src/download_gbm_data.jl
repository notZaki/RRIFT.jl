const gbm_study_uids = [
    "1.3.6.1.4.1.14519.5.2.1.4591.4001.278082550121070125285213632206",
    "1.3.6.1.4.1.14519.5.2.1.4591.4001.335353575986269052491315637674",
    "1.3.6.1.4.1.14519.5.2.1.4591.4001.365805576275232517344053939830",
    "1.3.6.1.4.1.14519.5.2.1.4591.4001.100057969162276274933613772317",
    "1.3.6.1.4.1.14519.5.2.1.4591.4001.269887096484012292940330991126",
    "1.3.6.1.4.1.14519.5.2.1.4591.4001.961791689281776173751323306588",
    "1.3.6.1.4.1.14519.5.2.1.4591.4001.763554173270318063812534542847",
    "1.3.6.1.4.1.14519.5.2.1.4591.4001.304604545029494418165835320551"
    ]

const alt_gbm_name = Dict(
    "TCGA-06-0185-1" => "1.3.6.1.4.1.14519.5.2.1.4591.4001.278082550121070125285213632206",
    "TCGA-06-0185-2" => "1.3.6.1.4.1.14519.5.2.1.4591.4001.335353575986269052491315637674",
    "TCGA-06-0185-3" => "1.3.6.1.4.1.14519.5.2.1.4591.4001.365805576275232517344053939830",
    "TCGA-06-0881-1" => "1.3.6.1.4.1.14519.5.2.1.4591.4001.100057969162276274933613772317",
    "TCGA-06-0881-2" => "1.3.6.1.4.1.14519.5.2.1.4591.4001.269887096484012292940330991126",
    "TCGA-06-1802-1" => "1.3.6.1.4.1.14519.5.2.1.4591.4001.961791689281776173751323306588",
    "TCGA-06-2570-1" => "1.3.6.1.4.1.14519.5.2.1.4591.4001.763554173270318063812534542847",
    "TCGA-06-5417-1" => "1.3.6.1.4.1.14519.5.2.1.4591.4001.304604545029494418165835320551"
    )

function download_invivo_studies(studies::AbstractVector{String}=gbm_study_uids; destination::AbstractString, overwrite = false)
    make_folder(destination)

    number_of_studies = length(studies)
    vfa_folders, dce_folders = [fill("", number_of_studies) for _=1:2]
    for (index, study) in enumerate(studies)
        (vfa_folders[index], dce_folders[index]) = download_vfa_and_dce(study, destination; overwrite=overwrite)
    end
    return (vfa_folders = vfa_folders, dce_folders = dce_folders)
end

download_invivo_studies(study::AbstractString; kwargs...) = download_invivo_studies([study]; kwargs...)

function download_vfa_and_dce(study_id, destination; overwrite = false)
    vfa_folder = joinpath(destination, study_id, "vfa")
    dce_folder = joinpath(destination, study_id, "dce")

    if isdir(vfa_folder) && isdir(dce_folder)
        if !overwrite
            return (vfa_folder, dce_folder)
        end
    end

    gbm_series = series(study = study_id)
    vfa_series = find_vfa_series(gbm_series)
    dce_series = find_dce_series(gbm_series)

    if !isdir(vfa_folder) || overwrite == true
        download_series(vfa_series, destination = vfa_folder)
    end

    if !isdir(dce_folder) || overwrite == true
        download_series(dce_series, destination = dce_folder)
    end
    return (vfa_folder, dce_folder)
end


find_vfa_series(series_dataframe) = find_in_description("MAP", series_dataframe)
find_dce_series(series_dataframe) = find_in_description("DYN", series_dataframe)

function find_in_description(word_to_find::AbstractString, series_dataframe)
    descriptions = series_dataframe.SeriesDescription
    found_indices = findall(occursin.(word_to_find, descriptions))
    if length(found_indices) < 1
        errant_id = series_dataframe.PatientID[1]
        errant_study = series_dataframe.StudyInstanceUID[1]
        error("No single series in $errant_id found with $word_to_find in their description.
              The full study UID is $errant_study\n")
    elseif length(found_indices) > 1
        errant_id = series_dataframe.PatientID[1]
        @warn "Multiple series in $errant_id found containing $word_to_find in their description.
              An arbitrary series will be used among the found series"
        found_index = found_indices[end]
    else
        found_index = found_indices[1]
    end
    found_series = series_dataframe.SeriesInstanceUID[found_index]
    return found_series
end

function download_series(series_id::AbstractString; destination)
    make_folder(destination; remove_existing = true)
    zip_file = joinpath(destination, "downloaded.zip")
    images(series = series_id, file = zip_file)
    unzip_command = `unzip -o $zip_file -d $destination`
    run(unzip_command)
    rm(zip_file)
    return destination
end

function download_invivo_masks(; destination, overwrite = false)
    if isdir(destination) && overwrite == false
        return destination
    end
    make_folder(destination; remove_existing = true)
    zip_file = joinpath(destination, "invivo_masks.zip")
    download("https://osf.io/uxe3p/download", zip_file)
    unzip_cmd = `unzip -o $zip_file -d $destination`
    run(unzip_cmd)
    rm(zip_file)
    rename_masks(folder = destination)
    return destination
end

function rename_masks(; folder)
    for (root, dirs, files) in walkdir(folder)
        for file in files
            if splitext(file)[2] == ".mat"
                oldfile = joinpath(root, file)
                oldfilename = splitext(file)[1]
                newfilename = alt_gbm_name[oldfilename] * ".mat"
                newfile = joinpath(root, newfilename)
                mv(oldfile, newfile)
            end
        end
    end
    return
end

function download_invivo_preprocessed(; destination, overwrite = false)
    if isdir(destination) && overwrite == false
        return destination
    end
    make_folder(destination; remove_existing = true)
    zip_file = joinpath(destination, "preprocessed.zip")
    download("https://osf.io/d3ext/download", zip_file)
    unzip_cmd = `unzip -o $zip_file -d $destination`
    run(unzip_cmd)
    rm(zip_file)
    return destination
end
