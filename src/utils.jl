function make_folder(desired_path; remove_existing = false)
    if !isdir(desired_path)
        mkpath(desired_path)
    elseif remove_existing
        rm(desired_path; recursive = true)
        mkpath(desired_path)
    end
    return desired_path
end

function ccc(x::AbstractVector, y::AbstractVector; lim=(0, Inf))
    @assert length(x) == length(y)
    mask = trues(length(x))
    for i in eachindex(x)
        mask[i] = !isnan(x[i]) & !isnan(y[i]) & (x[i]>lim[1]) & (x[i]<lim[2]) & (y[i]>lim[1]) & (y[i]<lim[2])
    end
    x = x[mask]
    y = y[mask]

    mx = mean(x)
    my = mean(y)
    sx = var(x) * (length(x) - 1) / length(x)
    sy = var(y) * (length(y) - 1) / length(y)
    sxy = sum((x .- mx) .* (y .- my)) / length(x)
    ccc_value = 2 * sxy / (sx + sy + (mx - mx)^2)
    return ccc_value
end

function apply_mask(; data, mask)
    @assert length(size(data)) >= length(size(mask))
    mask_indices = findall(mask)
    return data[mask_indices, :]
end

function negatives_to_zero!(x::AbstractArray)
    for i in eachindex(x)
        x[i] = x[i] < 0 ? 0 : x[i]
    end
    return x
end

function autocrop(data::AbstractArray; mask)
    xlim = findall(vec(sum(mask, dims = [2,3])) .> 0)
    ylim = findall(vec(sum(mask, dims = [1,3])) .> 0)
    return data[xlim, ylim, :]
end
function autocrop(data::NamedTuple; mask)
    for key in keys(data)
        data[key] = autocrop(data[key]; mask = mask)
    end
    return data
end
