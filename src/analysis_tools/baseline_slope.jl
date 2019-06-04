function linear_regression(x::Vector{<:Real}, y::Vector{<:Real})::Tuple # Substitutes linear fit --> much faster
    T=Float64
    x_mean::T = mean(x)
    y_mean::T = mean(y)

    num::T = 0.0
    nom::T = 0.0

    @fastmath for i in eachindex(x)
        x_res = (x[i] - x_mean)
        num += x_res * (y[i] - y_mean)
        nom += x_res*x_res
    end

    @fastmath slope::T = num/nom
    @fastmath offset::T = y_mean - slope * x_mean
    return offset, slope
end

function rms(v::Vector{<:Real})
    return sqrt(sum(v.^2)/length(v))
end



function determine_baseline_information(m::Measurement)
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    T = Float32
    sampling_rate = T(m.daq.sampling_rate)
    bl = Int(m.daq.baseline_length )
    bl_inv = T(1 / bl)
    # slope_length::Int = 100
    # slope_w1::UnitRange = 1:slope_length
    # slope_w2::UnitRange = bl - slope_length:bl
    # slope_length_inv::T = 1 / slope_length

    for f in inputfiles[1:1]
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
            d_baseline_rms = d_create(g_pd, "baseline_rms",    Float32, ((n_channel, n_events),(n_channel, n_events)), "chunk", chunksize)

            # n_events = 2000 # debugging

            # @fastmath @inbounds begin
            begin
                chunk_offsets = Array{T, 2}(undef, n_channel, chunksize[2])
                chunk_slopes = Array{T, 2}(undef, n_channel, chunksize[2])
                chunk_rms = Array{T, 2}(undef, n_channel, chunksize[2])
                chunk_pulses = Array{T, 3}(undef, n_samples, n_channel, chunksize[2])
                # @showprogress for evt_range in event_range_iterator(n_events, chunksize[2])
                for evt_range in event_range_iterator(n_events, chunksize[2])
                    @info evt_range
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
    return nothing
end

function flag_pileup_events(m, threshold = missing, channel::Int=1; overwrite = false)
    !exists(m,"Processed_data/baseline_slope") ? determine_baseline_information(m) : nothing
    n_total_events = get_number_of_events(m)
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    core::UInt8 = 1
    if overwrite
    for (fi, f) in enumerate(inputfiles)
        h5f = h5open(f, "r+")
        try
            g_pd  = g_open(h5f, "Processed_data")
            d_slopes = d_open(g_pd, "baseline_slope")
            T::Type = eltype(d_slopes)
            ismissing(threshold) ? d = T(-3*stdm(d_slopes[channel,:][1,:],0.0)) : d = T(threshold)
            n_channel, n_events = size(d_slopes)
            chunk_n_events = get_chunk(d_slopes)[end]
            if exists(g_pd, "pile_up_flags")
                o_delete(g_pd, "pile_up_flags")
            end
            d_pile_up_flag  = d_create(g_pd, "pile_up_flags", UInt8, ((n_events,),(n_events,)), "chunk", (chunk_n_events,) )

            evt_ranges = event_range_iterator(n_events, chunk_n_events)
            @fastmath @inbounds begin
                @showprogress for evt_range in evt_ranges
                    chunk_slopes::Array{T, 1} = d_slopes[channel, evt_range][1,:]
                    chunk_pileup::Array{UInt8, 1} = zeros(UInt8, length(evt_range))
                    for event in 1:length(evt_range)
                        chunk_pileup[event] = chunk_slopes[event] < d ? UInt8(1) : UInt8(0)
                    end
                    d_pile_up_flag[evt_range] = chunk_pileup
                end
            end

            close(h5f)
        catch err
            close(h5f)
            error(err)
        end
    end
    return nothing
end
# function get_baseline_slope_distribution(m::Measurement, n_events::Union{Int,Missing}=missing)::Vector{Float32}
#     T=Float32
#     ismissing(n_events) ? n_events  = get_number_of_events(m) : nothing
#     println("n_events: $n_events")
# 	inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
#     h5f = h5open(inputfiles[1], "r")
#         g_daq = g_open(h5f, "DAQ_Data")
#         d_daq_pulses = d_open(g_daq, "daq_pulses")
#         samplesize::T = 1e9 / m.daq.sampling_rate
#         slopes = zeros(T,n_events)
#         x = [ (i-1) * samplesize for i in 1:m.daq.baseline_length ]
#         daq_core_baselines =  d_daq_pulses[1:m.daq.baseline_length,1,1:n_events][:,1,:]
#         println(typeof(daq_core_baselines))
#     close(h5f)
#     @showprogress for ievt in 1:n_events
#         # slopes[ievt]= get_baseline_slope( d_daq_pulses[:,1,ievt][1:m.daq.baseline_length], samplesize )[2]
#         # slopes[ievt] = linear_regression(x, d_daq_pulses[:,1,ievt][1:m.daq.baseline_length])[2]
#         slopes[ievt] = linear_regression(x, daq_core_baselines[:,ievt])[2]
#     end
#     return slopes
# end
function get_baseline_slope(baseline::Vector, samplesize::Float32=4)::Tuple
    return linear_regression([ (i-1) * samplesize for i in eachindex(baseline) ], baseline)
end
