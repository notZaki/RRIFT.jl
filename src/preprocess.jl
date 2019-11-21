function load_preprocessed_mat(file)
    data = matread(file)
    return (t = data["t"], ct = data["ct"], crr = data["crr"], cp = data["cp"], relaxation = data["relaxation"], masks = data["masks"])
end

function preprocess_dicom_to_mat(; destination, dicom_folders, mask_folder, r1=3.3/1000, num_baseline=3, overwrite=false)
    @extract (vfa_folders, dce_folders) dicom_folders
    @assert length(vfa_folders) == length(dce_folders)
    make_folder(destination)
    num_folders = length(vfa_folders)
    matfiles = fill("", num_folders)
    for i = 1:num_folders
        study_id = split(vfa_folders[i], "/")[end-1]
        file = joinpath(destination, study_id * ".mat")
        if isfile(file) && overwrite == false
            matfiles[i] = file
            continue
        end
        computed = compute_concentration(vfa_folder = vfa_folders[i], dce_folder = dce_folders[i], r1 = r1, num_baseline = num_baseline)
        masks = get_mask(study = study_id, mask_folder = mask_folder)
        t = computed.t
        ct = apply_mask(data = computed.ct, mask = masks.tumour)
        cp = apply_mask(data = computed.ct, mask = masks.aif) ./ (1-0.4)
        crr = apply_mask(data = computed.ct, mask = masks.muscle)

        cp = mean(cp, dims=1)
        crr = mean(crr, dims=1)

        T1_tumour = apply_mask(data = computed.T1, mask = masks.tumour)
        T1_aif = apply_mask(data = computed.T1, mask = masks.aif)
        T1_muscle = apply_mask(data = computed.T1, mask = masks.muscle)
        M0_tumour = apply_mask(data = computed.M0, mask = masks.tumour)
        M0_aif = apply_mask(data = computed.M0, mask = masks.aif)
        M0_muscle = apply_mask(data = computed.M0, mask = masks.muscle)
        relaxation = Dict(
        "T1" => Dict("tumour" => T1_tumour, "muscle" => T1_muscle, "aif" => T1_aif),
        "M0" => Dict("tumour" => M0_tumour, "muscle" => M0_muscle, "aif" => M0_aif),
        )
        matfiles[i] = save_concentration_as_mat(file = file, t = t, ct = ct, crr = vec(crr), cp = vec(cp), relaxation = relaxation, masks = masks)
    end
    return matfiles
end

function get_mask(; study, mask_folder)
    aif_dir = joinpath(mask_folder, "AIF")
    muscle_dir = joinpath(mask_folder, "Muscle")
    tumour_dir = joinpath(mask_folder, "Tumour")
    filename = study * ".mat"

    aif = convert(BitArray{3}, matread(joinpath(aif_dir, filename))["mask"])
    muscle = convert(BitArray{3}, matread(joinpath(muscle_dir, filename))["mask"])
    tumour = convert(BitArray{3}, matread(joinpath(tumour_dir, filename))["mask"])

    return (aif = aif, muscle = muscle, tumour = tumour)
end

function save_concentration_as_mat(; file, t, ct, crr, cp, relaxation, masks)
    output_data = Dict(
    "t" => t,
    "ct" => ct,
    "crr" => crr,
    "cp" => cp,
    "relaxation" => relaxation,
    "masks" => masks
    )
    matwrite(file, output_data; compress = true)
    return file
end

function compute_concentration(; vfa_folder, dce_folder, r1=3.3/1000, num_baseline::Int=3)
    vfa = load_vfa_dicom(folder = vfa_folder)
    dce = load_dce_dicom(folder = dce_folder)

    relaxation_maps = fit_relaxation(:despot; vfa...).estimates
    concentration = signal_to_concentration(dce.signal, R10=1.0./relaxation_maps.T1, angle=dce.angle, TR=dce.TR, r1=r1, BAF=num_baseline)
    return (ct = concentration, t = dce.timepoints, relaxation_maps...)
end

function load_vfa_dicom(; folder)
    dicom_data = load_dicom_from_folder(folder)
    number_of_images = length(dicom_data)
    unique_flip_angles = unique(lookup.(dicom_data, "Flip Angle"))
    number_of_flip_angles = length(unique_flip_angles)
    number_of_slices = Int(number_of_images / number_of_flip_angles)

    dummy_image = lookup(dicom_data[1], "Pixel Data")'
    image_size = size(dummy_image)
    signal_data = zeros(image_size..., number_of_images)
    flip_angles = zeros(number_of_images)
    for dicom in dicom_data
        instance = lookup(dicom, "Instance Number")
        signal_data[:,:,instance] = lookup(dicom, "Pixel Data")'
        flip_angles[instance] = lookup(dicom, "Flip Angle")
    end
    signal_data = reshape(signal_data, (image_size..., number_of_slices, number_of_flip_angles))
    flip_angles = reshape(flip_angles, (number_of_slices, number_of_flip_angles))[1,:]
    @. flip_angles = deg2rad(flip_angles)
    TR = lookup(dicom_data[1], "Repetition Time")
    return (signal = signal_data, angles = flip_angles, TR = TR)
end

function load_dce_dicom(; folder, num_slices::Int = 16)
    dicom_data = load_dicom_from_folder(folder)
    number_of_images = length(dicom_data)
    number_of_timepoints = Int(number_of_images / num_slices)

    dummy_image = lookup(dicom_data[1], "Pixel Data")'
    image_size = size(dummy_image)
    signal_data = zeros(image_size..., number_of_images)
    timepoints = zeros(number_of_images)
    for dicom in dicom_data
        instance = lookup(dicom, "Instance Number")
        signal_data[:,:,instance] = lookup(dicom, "Pixel Data")'
        timepoints[instance] = lookup(dicom, "Trigger Time")
    end
    signal_data = reshape(signal_data, (image_size..., num_slices, number_of_timepoints))
    timepoints = reshape(timepoints, (num_slices, number_of_timepoints))[1,:] ./ 1000 ./ 60
    TR = lookup(dicom_data[1], "Repetition Time")
    flip_angle = deg2rad(lookup(dicom_data[1], "Flip Angle"))
    return (signal = signal_data, timepoints = vec(timepoints), TR = TR, angle = flip_angle)
end

function load_dicom_from_folder(folder)
    dicom_files = get_dicom_files_in_folder(folder = folder)
    number_of_files = length(dicom_files)
    dummy_data = dcm_parse(dicom_files[1])
    dicom_data = Vector{typeof(dummy_data)}(undef, number_of_files)
    for (index, file) in enumerate(dicom_files)
        dicom_data[index] = dcm_parse(file)
    end
    return dicom_data
end

function get_dicom_files_in_folder(; folder)
    files = readdir(folder)
    dicom_files = joinpath.(folder, files[is_dicom.(files)])
    return dicom_files
end

is_dicom(file::AbstractString) = splitext(file)[2] == ".dcm"
