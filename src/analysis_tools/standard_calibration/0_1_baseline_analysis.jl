@fastmath function linear_regression(x::Vector{<:Real}, y::Vector{<:Real})::Tuple # Substitutes linear fit --> much faster
    @assert length(x) == length(y) "x and y must have the same length."
    T=Float64
    x_mean::T = mean(x)
    y_mean::T = mean(y)
    num::T = 0.0
    nom::T = 0.0
    for i in eachindex(x)
        x_res = (x[i] - x_mean)
        num += x_res * (y[i] - y_mean)
        nom += x_res*x_res
    end
    slope::T = num/nom
    offset::T = y_mean - slope * x_mean
    return offset, slope
end
@fastmath function linear_regression(x::Vector{T}, y::Vector{T})::Tuple{T, T} where {T <: AbstractFloat} # Substitutes linear fit --> much faster
    @assert length(x) == length(y) "x and y must have the same length."
    x_mean::T = mean(x)
    y_mean::T = mean(y)
    num::T = 0.0
    nom::T = 0.0
    @inbounds for i in eachindex(x)
        x_res = (x[i] - x_mean)
        num += x_res * (y[i] - y_mean)
        nom += x_res * x_res
    end
    slope::T = num / nom
    offset::T = y_mean - slope * x_mean
    return offset, slope
end

function rms(v::Vector{<:Real})
    return sqrt(sum(v.^2)/length(v))
end

function determine_baseline_information(m::Measurement; n_sigma::Real = 1.0)
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    determine_baseline_information(inputfiles, m, n_sigma = n_sigma)
end
function determine_baseline_information(filename::String, m::Measurement; n_sigma::Real = 1.0)
    determine_baseline_information([filename], m,  n_sigma = n_sigma)
end
function determine_baseline_information(filenames::Vector{String}, m::Measurement; n_sigma::Real = 1.0)
    inputfiles = filenames
    T = Float32
    sampling_rate = T(m.daq.sampling_rate)
    bl = Int(m.daq.baseline_length )
    bl_inv = T(1 / bl)

    for f in inputfiles
        new_pulse_format = is_new_pulse_format(f)
        if !new_pulse_format error("Measurement $(m.name) does not have new pulse format. Old one ist not implemented yet.") end
        h5f = h5open(f, "r+")
        try
            g_daq = g_open(h5f, "DAQ_Data")
            d_daq_pulses = d_open(g_daq, "daq_pulses")
            daq_pulse_T = eltype(d_daq_pulses)
            chunksize_daq = get_chunk(d_daq_pulses)
            n_samples, n_channel, n_events = new_pulse_format ? size(d_daq_pulses) : (size(d_daq_pulses, 2), size(d_daq_pulses, 1), size(d_daq_pulses, 3))
            chunksize = n_channel, chunksize_daq[3]
            g_pd = exists(h5f, "Processed_data") ? g_open(h5f, "Processed_data") : g_create(h5f, "Processed_data")
            if exists(g_pd, "baseline_offset") o_delete(g_pd, "baseline_offset") end
            if exists(g_pd, "baseline_slope") o_delete(g_pd, "baseline_slope") end
            if exists(g_pd, "baseline_rms") o_delete(g_pd, "baseline_rms") end
            d_baseline_offset = d_create(g_pd, "baseline_offset", Float32, ((n_channel, n_events),(n_channel, n_events)), "chunk", chunksize)
            d_baseline_slope  = d_create(g_pd, "baseline_slope",  Float32, ((n_channel, n_events),(n_channel, n_events)), "chunk", chunksize)
            d_baseline_rms    = d_create(g_pd, "baseline_rms",    Float32, ((n_channel, n_events),(n_channel, n_events)), "chunk", chunksize)

            # n_events = 2000 # debugging

            # @fastmath @inbounds begin
            begin
                chunk_offsets = Array{T, 2}(undef, n_channel, chunksize[2])
                chunk_slopes = Array{T, 2}(undef, n_channel, chunksize[2])
                chunk_rms = Array{T, 2}(undef, n_channel, chunksize[2])
                chunk_pulses = Array{T, 3}(undef, n_samples, n_channel, chunksize[2])
                @showprogress for evt_range in event_range_iterator(n_events, chunksize[2])
                # for evt_range in event_range_iterator(n_events, chunksize[2])
                    # @info evt_range
                    levtr::Int = length(evt_range)
                    chunk_pulses[:, :, 1:levtr]  = d_daq_pulses[:, :, evt_range]
                    for i in 1:levtr
                        for ichn in 1:n_channel
                            baseline::Vector{T} = T.( chunk_pulses[1:bl, ichn, i] )
                            # slope::T = slope_length_inv * (mean(baseline[slope_w2]) - mean(baseline[slope_w1]))
                            slope = linear_regression(collect(0:bl-1),baseline)[2]
                            μ::T = mean(baseline)
                            σ::T = std(baseline.-μ)
                            chunk_slopes[ichn, i] = slope
                            chunk_offsets[ichn, i] = μ
                            chunk_rms[ichn, i] = σ
                        end
                    end
                    d_baseline_offset[:, evt_range] = chunk_offsets[:, 1:levtr]
                    d_baseline_slope[:, evt_range]  = chunk_slopes[:, 1:levtr]
                    d_baseline_rms[:, evt_range]    = chunk_rms[:, 1:levtr]
                end
            end

            close(h5f)
        catch err
            close(h5f)
            error(err)
        end
    end
    baseline_quality_plots(m, n_sigma = n_sigma)
    return nothing
end

function unflag_events(m)
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    for (fi, f) in enumerate(inputfiles)
        h5f = h5open(f, "r+")
        try
            g_pd  = g_open(h5f, "Processed_data")
            d_pile_up_flag = if exists(g_pd, "event_flags")
                d_open(g_pd, "event_flags")
            else
                d_create(g_pd, "event_flags", UInt8, ((n_events,),(n_events,)), "chunk", (chunk_n_events,) )
            end

            d_pile_up_flag[:] = HealthyEvent

            close(h5f)
        catch err
            close(h5f)
            error(err)
        end
    end
    return nothing
end

function flag_pileup_events(m, channel::Int = 1, n_sigma::Real = 1)
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    n_pile_up_events::Int = 0
    n_total_events::Int = get_number_of_events(m)
    for (fi, f) in enumerate(inputfiles)
        !exists(f, "Processed_data/baseline_slope") ? determine_baseline_information(m) : nothing
        h5f = h5open(f, "r+")
        try
            g_pd  = g_open(h5f, "Processed_data")
            d_slopes = d_open(g_pd, "baseline_slope")
            T::DataType = eltype(d_slopes)
            σ::T = stdm(d_slopes[channel,:][1,:],0.0)
            d::T = T(n_sigma * σ)
            n_channel::Int, n_events::Int = size(d_slopes)
            chunk_n_events::Int = get_chunk(d_slopes)[end]
            d_pile_up_flag_exists::Bool = exists(g_pd, "event_flags")
            d_pile_up_flag = if d_pile_up_flag_exists
                d_open(g_pd, "event_flags")
            else
                d_create(g_pd, "event_flags", UInt8, ((n_events,),(n_events,)), "chunk", (chunk_n_events,) )
            end
            if !d_pile_up_flag_exists d_pile_up_flag[:] .= UInt8(0) end

            slopes::Vector{T} = read(d_slopes)[channel, :]
            pile_up_flags::Vector{UInt8} = d_pile_up_flag[:]

            @inbounds for i in 1:n_events
                if slopes[i] < -d
                    if pile_up_flags[i] & PileUpEvent == 0
                        pile_up_flags[i] += PileUpEvent::UInt8 #  = UInt8(1)
                    end
                end
            end
            d_pile_up_flag[:] = pile_up_flags

            n_pile_up_events += length(findall(x -> (x & PileUpEvent) == PileUpEvent, pile_up_flags))

            close(h5f)
        catch err
            println(err)
            close(h5f)
            error(err)
        end
    end

    @info "$(m.name) - $n_pile_up_events pile up events ($(round(100 * n_pile_up_events / n_total_events, digits = 3)) %)"
    return nothing
end
function flag_rising_baseline_events(m, channel::Int = 1, n_sigma::Real = 1)
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    n_pile_up_events::Int = 0
    n_total_events::Int = get_number_of_events(m)
    for (fi, f) in enumerate(inputfiles)
        !exists(f, "Processed_data/baseline_slope") ? determine_baseline_information(m) : nothing
        h5f = h5open(f, "r+")
        try
            g_pd  = g_open(h5f, "Processed_data")
            d_slopes = d_open(g_pd, "baseline_slope")
            T::DataType = eltype(d_slopes)
            σ::T = stdm(d_slopes[channel,:][1,:],0.0)
            d::T = T(n_sigma * σ)
            n_channel::Int, n_events::Int = size(d_slopes)
            chunk_n_events::Int = get_chunk(d_slopes)[end]
            d_pile_up_flag_exists::Bool = exists(g_pd, "event_flags")
            d_pile_up_flag = if d_pile_up_flag_exists
                d_open(g_pd, "event_flags")
            else
                d_create(g_pd, "event_flags", UInt8, ((n_events,),(n_events,)), "chunk", (chunk_n_events,) )
            end
            if !d_pile_up_flag_exists d_pile_up_flag[:] .= UInt8(0) end

            slopes::Vector{T} = read(d_slopes)[channel, :]
            pile_up_flags::Vector{UInt8} = d_pile_up_flag[:]

            @inbounds for i in 1:n_events
                if slopes[i] > d
                    if pile_up_flags[i] & RisingBaselineEvent == 0
                        pile_up_flags[i] += RisingBaselineEvent::UInt8 #  = UInt8(1)
                    end
                end
            end
            d_pile_up_flag[:] = pile_up_flags

            n_pile_up_events += length(findall(x -> (x & RisingBaselineEvent) == RisingBaselineEvent, pile_up_flags))

            close(h5f)
        catch err
            close(h5f)
            error(err)
        end
    end

    @info "$(m.name) - $n_pile_up_events pile up events ($(round(100 * n_pile_up_events / n_total_events, digits = 3)) %)"
    return nothing
end

function get_baseline_slope(baseline::Vector, samplesize::Float32=4)::Tuple
    return linear_regression([ (i-1) * samplesize for i in eachindex(baseline) ], baseline)
end

function baseline_quality_plots(m::Measurement; n_sigma::Real = 1.0)
    files = gather_absolute_paths_to_hdf5_input_files(m)
    T::DataType = get_eltype_of_dataset(m, "Processed_data", "baseline_offset")
    n_events::Int = get_number_of_events(m)
    n_channel::Int = get_number_of_channel(m)
    # @info "\t $n_events Events"

    baseline_offsets::Array{T, 2}   = Array{T, 2}(undef, n_events, n_channel)
    baseline_rms::Array{T, 2}       = Array{T, 2}(undef, n_events, n_channel)
    baseline_slopes::Array{T, 2}    = Array{T, 2}(undef, n_events, n_channel)

    last_event_idx::Int = 0
    for fn in files
        h5open(fn, "r") do h5f
            g_pd = g_open(h5f, "Processed_data")
            d_blμ = d_open(g_pd, "baseline_offset")
            d_blσ = d_open(g_pd, "baseline_rms")
            d_bls = d_open(g_pd, "baseline_slope")
            n_events_in_file::Int = size(d_blμ, 2)
            evt_range::UnitRange{Int} = last_event_idx+1:last_event_idx+n_events_in_file
            baseline_offsets[evt_range, :]  = read(d_blμ)'
            baseline_rms[evt_range, :]      = read(d_blσ)'
            baseline_slopes[evt_range, :]   = read(d_bls)'
            last_event_idx += n_events_in_file
        end
    end

    for ichn in 1:n_channel
        begin # baseline offset
            bl_μ::T, bl_σ::T = mean_and_std(baseline_offsets[:, ichn])
            h_offsets = Histogram(0.0:1:2^14)
            append!(h_offsets,  baseline_offsets[:, ichn])
            # h_offsets = fit(Histogram, baseline_offsets[:, ichn], bl_μ-10bl_σ:bl_σ/100:bl_μ+10bl_σ, closed=:left)
            p_h_bl_offsets = plot(h_offsets, st=:step, yscale = :log10, label = "", xlabel = "Baseline mean / ADC counts", ylabel = "#Events")
            vline!([bl_μ - n_sigma * bl_σ, bl_μ + n_sigma * bl_σ], lw=2.0, label = "Mean ± $(n_sigma)σ", lc = :red)
        end

        begin # baseline rms
            bl_rms_μ::T, bl_rms_σ::T = mean_and_std(baseline_rms[:, ichn])
            h_bl_rms = fit(Histogram, baseline_rms[:, ichn], 0:bl_rms_σ/500:bl_rms_μ + 2bl_rms_σ)
            p_h_bl_rms = plot(h_bl_rms, st=:step, yscale = :log10, xlabel = "Basline RMS / ADC counts", ylabel = "#Events", label = "")
            vline!([bl_rms_μ], label = "Mean", lw = 2.0, lc = :red)
            xlims!(0, bl_rms_μ + 1bl_rms_σ)
        end

        begin # baseline slope
            bl_slope_μ::T, bl_slope_σ::T = mean_and_std(baseline_slopes[:, ichn])
            h_bl_slopes = fit(Histogram, baseline_slopes[:, ichn], -25bl_slope_σ:bl_slope_σ/100:25bl_slope_σ)
            p_h_bl_slopes = plot(h_bl_slopes, st=:step, yscale = :log10, xlabel = "Basline slopes / (ADC counts / sample)", ylabel = "#Events", label = "")
            vline!([bl_slope_μ - n_sigma * bl_slope_σ, bl_slope_μ + n_sigma * bl_slope_σ], label = "Mean ± $(n_sigma)σ", lw = 2.0, lc = :red)
        end

        p_summary = plot(
            p_h_bl_offsets, p_h_bl_rms, p_h_bl_slopes,
            layout = (@layout [a b c]), size = (1600, 500),
            title = "Channel $ichn",
            legendfontsize = 14, guidefontsize = 14, tickfontsize = 12, titlefontsize = 14
        )

        savefig(m, p_summary, "0_1_baseline_quality", "0_1_baseline_quality_chn$ichn", fmt=:png )
    end
    nothing
end
