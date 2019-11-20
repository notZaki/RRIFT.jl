function fit_errm(; t::AbstractVector, crr::AbstractVector, ct::AbstractArray, mask=true)
    @assert length(t) == length(crr) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    x1, x2, x3, x4 = (fill(NaN, volume_size...) for _=1:4)
    resolved_mask = resolve_mask_size(mask, volume_size)

    M = zeros(num_timepoints, 4)
    M[:,1] = cumul_integrate(t, crr, TrapezoidalFast())
    M[:,2] = cumul_integrate(t, M[:,1], TrapezoidalFast())
    M[:,4] = crr
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        cur_ct = cumul_integrate(t, ct[idx, :], TrapezoidalFast())
        M[:, 3] = -cumul_integrate(t, cur_ct, TrapezoidalFast())
        x1[idx], x2[idx], x3[idx], x4[idx] = M \ cur_ct
    end
    # Apply correction because fit actually returns: [mess, mess, kep, vp/kt_rr]
    a = x1 ./ x4
    b = x2 ./ x4
    root_term = @. a^2 - 4*b
    @. root_term[root_term < 0] = 0
    kep_rr = @. (a - sqrt(root_term)) / 2
    rel_kt = @. x4 * (a - x3 - kep_rr)
    rel_ve = @. rel_kt * kep_rr / x3
    kep = x3
    rel_vp = x4
    return (rel_kt = rel_kt, rel_ve = rel_ve, rel_vp = rel_vp, kep = kep, kep_rr = kep_rr)
end

function fit_cerrm(; t::AbstractVector, crr::AbstractVector, ct::AbstractArray, mask=true, kep_rr=0.0)
    @assert length(t) == length(crr) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    x1, rel_vp, kep = (fill(NaN, volume_size...) for _=1:3)
    resolved_mask = resolve_mask_size(mask, volume_size)

    if kep_rr <= 0
        rrm_estimates = fit_errm(t=t, crr=crr, ct=ct, mask=mask)
        viable_estimates = positive_only_mask(rrm_estimates)
        kep_rr = interquartile_mean(rrm_estimates.kep_rr[viable_estimates])
    end

    M = zeros(num_timepoints, 3)
    crr_int = cumul_integrate(t, crr, TrapezoidalFast())
    crr_int2 = cumul_integrate(t, crr_int, TrapezoidalFast())
    M[:,1] = crr_int + kep_rr * crr_int2
    M[:,2] = crr + kep_rr * crr_int
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        cur_ct = cumul_integrate(t, ct[idx, :], TrapezoidalFast())
        M[:, 3] = -cumul_integrate(t, cur_ct, TrapezoidalFast())
        x1[idx], rel_vp[idx], kep[idx] = M \ cur_ct
    end
    rel_kt = @. x1 - kep * rel_vp
    rel_ve = @. rel_kt * kep_rr / kep
    return (rel_kt = rel_kt, rel_ve = rel_ve, rel_vp = rel_vp, kep = kep, kep_rr = kep_rr)
end


