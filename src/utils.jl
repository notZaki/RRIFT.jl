function make_folder(desired_path; remove_existing = false)
    if !isdir(desired_path)
        mkpath(desired_path)
    elseif remove_existing
        rm(desired_path; recursive = true)
        mkpath(desired_path)
    end
    return desired_path
end

function apply_mask(; data, mask)
    @assert length(size(data)) >= length(size(mask))
    mask_indices = findall(mask)
    return data[mask_indices, :]
end

